#!/usr/bin/env python
"""Convert one output type (e.g. genotype/) of per-sim .tsv files into a Parquet dataset.

Sim parameters are parsed from the file names and added as columns. Wide per-individual
tables (generation + one column per individual/group) are stored long: one row per
generation x ind (or group) with a single value column named after the type.
The dataset is partitioned into folders by the parameters of the --layout (constant env:
shift_size / reg / n; moving optimum: shift_speed / reg / n; fluctuating env: fluc_opt /
fluc_interval / reg / n), so readers can
skip whole folders; each folder holds part files of up to --sims-per-part sims:

    <out_dir>/<type>/shift_size=9.4/reg=ave/n=10/part-0.parquet

A _SUCCESS.json summary is written last, only if every input file is accounted for.

Usage:
    python convert_parquet.py --raw-dir DATASET_DIR --type genotype --out-dir PARQUET_DIR

slurm_job.sh runs it per job on that job's sims only, with --part-prefix and --summary so
files from different jobs land side by side in the same partition folders.
"""

import argparse
import contextlib
import json
import multiprocessing as mp
import os
import re
import shutil
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault("POLARS_MAX_THREADS", "1")  # parallelism comes from the worker processes
import polars as pl  # noqa: E402

# file-name layouts: the named groups are the sim parameters (added as columns or used as
# partition folders); partition columns are stored in folder names, not inside the files
LAYOUTS = {
    # constant environment: n_2_psi_0.0_sz_9.4_reg_ave_<rep>_<ID>.tsv
    "cons": (re.compile(r"^n_(?P<n>\d+)_psi_(?P<psi>-?[0-9.]+)_sz_(?P<shift_size>-?[0-9.]+)"
                        r"_reg_(?P<reg>fit|ext|ave)_(?P<rep>\d+)_(?P<ID>\d+)\.tsv$"),
             ["shift_size", "reg", "n"]),
    # constant environment, optimum moving at shift_speed per tick: n_2_psi_0.0_sp_0.01_reg_ave_<rep>_<ID>.tsv
    "moving": (re.compile(r"^n_(?P<n>\d+)_psi_(?P<psi>-?[0-9.]+)_sp_(?P<shift_speed>-?[0-9.]+)"
                          r"_reg_(?P<reg>fit|ext|ave)_(?P<rep>\d+)_(?P<ID>\d+)\.tsv$"),
               ["shift_speed", "reg", "n"]),
    # fluctuating environment: n_2_psi_0.0_sz_20_reg_ave_gap_20_<rep>_<ID>.tsv
    "fluc": (re.compile(r"^n_(?P<n>\d+)_psi_(?P<psi>-?[0-9.]+)_sz_(?P<fluc_opt>-?[0-9.]+)"
                        r"_reg_(?P<reg>fit|ext|ave)_gap_(?P<fluc_interval>\d+)_(?P<rep>\d+)_(?P<ID>\d+)\.tsv$"),
             ["fluc_opt", "fluc_interval", "reg", "n"]),
}

# wide types -> (name of the id column, name of the value column)
LONG_TYPES = {
    "phenotype": ("ind", "phenotype"),
    "genotype": ("ind", "genotype"),
    "fitness": ("ind", "fitness"),
    "copied_val": ("ind", "copied_val"),
    "group_tag": ("ind", "group_tag"),
    "vars": ("ind", "var"),  # each group's variance, repeated for every member
}

# SLiM columns that differ from the parameter columns only by case (metadata N vs n,
# allele id vs ID) are renamed so the data also works in case-insensitive tools (SQL)
RENAMES = {"N": "N_pop", "id": "mut_id"}

# column types of the non-long types (after RENAMES), as unified over the full 59,400-sim
# run. Each sim job converts its own sims, so types are fixed here rather than unified
# afterwards; this keeps part files from different jobs readable as one dataset.
F, I, S = pl.Float64, pl.Int64, pl.String
FIXED_SCHEMAS = {
    "data": {
        "tick": I, "optimum": F, "seg_site_pop": I, "num_fixed_mut": I, "z_opt": F, "z_opt_closest": F,
        "g_opt": F, "s": F, "b_gz": F, "dz1_pred": F, "dg_sel": F, "dg_obs": F,
        "Ne_heterozygosity": F, "H": F, "PI": F, "sfs": S, "z_mean": F, "z_sd": F,
        "g_mean": F, "g_sd": F, "fitness_mean": F, "fitness_sd": F, "parent_z_mean": F,
        "parent_z_sd": S, "parent_g_mean": F, "parent_g_sd": S, "group_vars_mean": F,
        "group_vars_sd": F, "num_seg_mean": F, "num_seg_sd": F, "abs_effect_mean": F,
        "abs_effect_sd": F, "mean_abs_effect_mean": F, "mean_abs_effect_sd": F,
        "mean_effect_mean": F, "mean_effect_sd": F,
    },
    "metadata": {
        "N_pop": I, "Vs": I, "Va0_sqrt": F, "optimum": F, "seed": I, "end_time": I,
        "shift_time": I, "maxL": I, "U": F, "stop_tick": I, "stopped_early": S,
        "stop_fit_ratio": F, "fit_ratio": F, "simulation_time": F,
    },
    "allele": {
        "gen": I, "mut_id": I, "position": I, "effect_size": F, "counts": I,
        "frequencies": F, "isFixed": S, "originTick": I,
    },
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--layout", choices=LAYOUTS, default="cons", help="file-name layout")
    p.add_argument("--raw-dir", required=True, help="dataset folder holding the <type>/ subfolders")
    p.add_argument("--type", required=True, nargs="+",
                   help="output type(s), e.g. genotype; several are converted in one process")
    p.add_argument("--out-dir", required=True, help="parquet root; writes <out-dir>/<type>/")
    p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    p.add_argument("--sims-per-part", type=int, default=100,
                   help="sims per parquet part file (memory per worker grows with this)")
    p.add_argument("--part-prefix", default="",
                   help="added to part file names (part-<prefix><k>.parquet) so files written "
                        "by different jobs into the same partition folder do not collide")
    p.add_argument("--summary-dir", help="write each summary to <dir>/<type>.json "
                                         "(default <out-dir>/<type>/_SUCCESS.json)")
    return p.parse_args()


def scan(type_dir, name_rx):
    """Return (sim records, unparseable names, empty files) for the .tsv files in type_dir."""
    sims, bad, empty = [], [], []
    with os.scandir(type_dir) as it:
        for e in it:
            if not e.name.endswith(".tsv"):
                continue
            m = name_rx.match(e.name)
            if not m:
                bad.append(e.name)
            elif e.stat().st_size == 0:
                empty.append(e.name)
            else:
                g = m.groupdict()  # other parameters stay strings, as in the folder names
                sims.append(g | {"path": e.path, "n": int(g["n"]), "psi": float(g["psi"]),
                                 "rep": int(g["rep"]), "ID": int(g["ID"])})
    return sims, bad, empty


def read_one(sim, sim_type, meta_cols):
    """Read one sim's tsv into a DataFrame with the parameter columns added."""
    if sim_type == "vars":
        # older SLiM runs wrote only N/n header columns for N values per row; the values
        # are per individual (index 0..N-1, as in fitness), so ignore the header
        df = pl.read_csv(sim["path"], separator="\t", has_header=False, skip_rows=1,
                         infer_schema_length=None)
        df.columns = ["generation"] + [str(i) for i in range(df.width - 1)]
    else:
        df = pl.read_csv(sim["path"], separator="\t", infer_schema_length=None)
    # blank lines (allele files have one after the header) come in as all-null rows
    df = df.filter(~pl.all_horizontal(pl.all().is_null()))
    if sim_type in LONG_TYPES:
        key, value = LONG_TYPES[sim_type]
        df = (df.unpivot(index="generation", variable_name=key, value_name=value)
                .drop_nulls(value)
                .with_columns(pl.col("generation").cast(pl.Int32),
                              pl.col(key).cast(pl.Int32),
                              pl.col(value).cast(pl.Float64)))
    else:
        # parameters from the file name win over SLiM's copies of them
        df = df.drop([c for c in df.columns if c in meta_cols])
        df = df.rename({c: RENAMES.get(c, c + "_slim") for c in df.columns
                        if c.lower() in {m.lower() for m in meta_cols}})
        schema = FIXED_SCHEMAS.get(sim_type, {})
        # columns not in the schema keep their inferred type; missing ones come in as null
        df = df.with_columns([pl.col(c).cast(dt) if c in df.columns else pl.lit(None, dt).alias(c)
                              for c, dt in schema.items()])
        df = df.select([*schema, *[c for c in df.columns if c not in schema]])
    return df.with_columns(pl.lit(sim["psi"]).cast(pl.Float64).alias("psi"),
                           pl.lit(sim["rep"]).cast(pl.Int32).alias("rep"),
                           pl.lit(sim["ID"]).cast(pl.Int32).alias("ID"))


def write_part(job):
    """Read a batch of sims from one partition and write them as one parquet file."""
    sims, sim_type, meta_cols, path = job
    frames, failed = [], []
    for s in sims:
        try:
            frames.append(read_one(s, sim_type, meta_cols))
        except Exception as e:  # noqa: BLE001 - report and keep going
            failed.append((os.path.basename(s["path"]), f"{type(e).__name__}: {e}"[:300]))
    if not frames:
        return path, [], failed, 0
    df = pl.concat(frames, how="diagonal_relaxed")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.write_parquet(path, compression="zstd", statistics=True)
    return path, df["ID"].unique().to_list(), failed, df.height


def unify_schemas(paths):
    """Rewrite part files whose column types differ, e.g. an all-integer column that is
    decimal in other parts, so every part shares one schema."""
    schemas = {p: pl.read_parquet_schema(p) for p in paths}
    target = {}
    for sch in schemas.values():
        for col, dt in sch.items():
            if col not in target:
                target[col] = dt
            elif target[col] != dt:
                if target[col].is_numeric() and dt.is_numeric():
                    target[col] = pl.Float64 if (target[col].is_float() or dt.is_float()) else pl.Int64
                else:
                    target[col] = pl.String
    n_rewritten = 0
    for p, sch in schemas.items():
        if dict(sch) != target:
            df = pl.read_parquet(p)
            df = df.with_columns([pl.lit(None).alias(c) for c in target if c not in df.columns])
            df = df.select([pl.col(c).cast(dt) for c, dt in target.items()])
            df.write_parquet(p, compression="zstd", statistics=True)
            n_rewritten += 1
    return target, n_rewritten


def n_data_lines(path):
    """Non-blank lines after the header."""
    with open(path) as f:
        return max(0, sum(1 for line in f if line.strip()) - 1)


class ConvertError(Exception):
    pass


def convert_type(a, sim_type):
    """Convert one type; raises ConvertError if any sim with data is not in the parquet."""
    t0 = time.time()
    type_dir = os.path.join(a.raw_dir, sim_type)
    out = os.path.join(a.out_dir, sim_type)

    name_rx, partition_cols = LAYOUTS[a.layout]
    meta_cols = list(name_rx.groupindex)
    sims, bad, empty = scan(type_dir, name_rx)
    print(f"[{sim_type}] {len(sims)} files to convert, {len(bad)} unrecognized names, {len(empty)} empty",
          flush=True)
    if not sims:
        raise ConvertError(f"No convertible .tsv files in {type_dir}")
    if os.path.exists(out):
        shutil.rmtree(out)

    # group by partition, keep rows ordered by psi / rep so readers can skip row groups
    groups = defaultdict(list)
    for s in sims:
        groups[tuple(s[c] for c in partition_cols)].append(s)
    jobs = []
    for key, members in sorted(groups.items()):
        members.sort(key=lambda s: (s["psi"], s["rep"], s["ID"]))
        folder = os.path.join(out, *[f"{c}={v}" for c, v in zip(partition_cols, key)])
        for k in range(0, len(members), a.sims_per_part):
            path = os.path.join(folder, f"part-{a.part_prefix}{k // a.sims_per_part}.parquet")
            jobs.append((members[k:k + a.sims_per_part], sim_type, meta_cols, path))

    found, failed, parts, n_rows = set(), [], [], 0
    # workers=1 runs in this process. Otherwise spawn, not fork: forking after polars has
    # started its threads (e.g. for an earlier type in this process) can deadlock the workers
    with (contextlib.nullcontext() if a.workers == 1 else
          ProcessPoolExecutor(max_workers=a.workers, mp_context=mp.get_context("spawn"))) as ex:
        results = (map(write_part, jobs) if ex is None else
                   (f.result() for f in as_completed([ex.submit(write_part, j) for j in jobs])))
        for i, (path, ids, fails, rows) in enumerate(results, 1):
            found.update(ids)
            failed.extend(fails)
            n_rows += rows
            if rows:
                parts.append(path)
            if i % max(1, len(jobs) // 10) == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} parts ({(time.time() - t0) / 60:.1f} min)", flush=True)

    if not parts:
        raise ConvertError(f"[{sim_type}] no parquet written; {len(failed)} read errors, e.g. {failed[:1]}")
    schema, n_rewritten = unify_schemas(parts)
    print(f"[{sim_type}] {len(parts)} part files, {n_rewritten} rewritten to a common schema", flush=True)

    # every sim with data rows must show up in the parquet
    on_disk = set(pl.scan_parquet(os.path.join(out, "**", "*.parquet"), hive_partitioning=True)
                  .select(pl.col("ID").unique()).collect()["ID"].to_list())
    missing = [s for s in sims if s["ID"] not in on_disk]
    header_only = [os.path.basename(s["path"]) for s in missing if n_data_lines(s["path"]) == 0]
    lost = [os.path.basename(s["path"]) for s in missing if n_data_lines(s["path"]) > 0]

    summary = {
        "type": sim_type,
        "raw_dir": type_dir,
        "files_converted": len(sims),
        "sims_in_parquet": len(on_disk),
        "rows": n_rows,
        "columns": {c: str(dt) for c, dt in schema.items()} | {c: "partition" for c in partition_cols},
        "header_only_files": header_only,
        "empty_files": empty,
        "unrecognized_names": bad,
        "read_errors": failed,
        "lost_files": lost,
        "minutes": round((time.time() - t0) / 60, 2),
    }
    print(json.dumps({k: (v if not isinstance(v, list) or len(v) < 10 else f"{len(v)} items")
                      for k, v in summary.items()}, indent=2))
    if lost:
        raise ConvertError(f"[{sim_type}] {len(lost)} files with data are missing from the parquet, e.g. {lost[0]}")
    summary_path = (os.path.join(a.summary_dir, f"{sim_type}.json") if a.summary_dir
                    else os.path.join(out, "_SUCCESS.json"))
    os.makedirs(os.path.dirname(os.path.abspath(summary_path)), exist_ok=True)
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[{sim_type}] done: {len(on_disk)} sims, {n_rows:,} rows")


def main():
    a = parse_args()
    failed = []
    for sim_type in a.type:
        try:
            convert_type(a, sim_type)
        except ConvertError as e:
            print(f"ERROR: {e}", file=sys.stderr, flush=True)
            failed.append(sim_type)
    if failed:
        sys.exit(f"conversion failed for: {' '.join(failed)}")


if __name__ == "__main__":
    main()
