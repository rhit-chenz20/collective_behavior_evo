#!/usr/bin/env python
"""Check that vars, copied_val and fitness can be recomputed from the other outputs, before
deleting them. Read-only. For a random sample of sims it recomputes, per generation:

    vars       = var(genotype) within the individual's group (genotype + group_tag)
    copied_val = the value the model copies, from the group's genotypes:
                 ave: mean of the other members; fit: member closest to the optimum;
                 ext: member farthest from the group mean. Ties (e.g. ext with n=2, where both
                 members are equally extreme and SLiM takes whichever came first in the random
                 shuffle) count as a match if the stored value is any of the tied members.
    fitness    = dnorm(phenotype, optimum, sqrt(Vs)), optimum = z_opt + z_mean from data

and reports how many stored values match (tolerance covers the tsv rounding).

Usage:
    python verify_derived.py --parquet-dir <dataset>/parquet [--sims 30]
"""

import argparse
import math
import os
import random

import polars as pl

KEY = ["ID", "generation", "ind"]


def scan(root, t, ids):
    return (pl.scan_parquet(os.path.join(root, t, "**", "*.parquet"), hive_partitioning=True)
            .filter(pl.col("ID").is_in(ids)))


def close(a, b, rtol=1e-3, atol=1e-6):
    return (a - b).abs() <= atol + rtol * b.abs()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--parquet-dir", required=True)
    p.add_argument("--sims", type=int, default=30, help="number of sims to check")
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args()
    root = a.parquet_dir

    all_ids = (pl.scan_parquet(os.path.join(root, "metadata", "**", "*.parquet"), hive_partitioning=True)
               .select("ID").collect()["ID"].to_list())
    ids = random.Random(a.seed).sample(all_ids, min(a.sims, len(all_ids)))

    g = scan(root, "genotype", ids).select(*KEY, "genotype", "reg", "n", "psi")
    t = scan(root, "group_tag", ids).select(*KEY, pl.col("group_tag").cast(pl.Int64))
    z = scan(root, "phenotype", ids).select(*KEY, "phenotype")
    stored = (scan(root, "vars", ids).select(*KEY, "var")
              .join(scan(root, "copied_val", ids).select(*KEY, "copied_val"), on=KEY)
              .join(scan(root, "fitness", ids).select(*KEY, "fitness"), on=KEY))
    opt = (scan(root, "data", ids)
           .select("ID", pl.col("tick").cast(pl.Int32).alias("generation"),
                   (pl.col("z_opt") + pl.col("z_mean")).alias("optimum")))
    vs = scan(root, "metadata", ids).select("ID", "Vs")

    df = (g.join(t, on=KEY).join(z, on=KEY).join(stored, on=KEY)
          .join(opt, on=["ID", "generation"]).join(vs, on="ID")
          .collect())
    grp = ["ID", "generation", "group_tag"]
    df = df.with_columns(
        var_pred=pl.col("genotype").var(ddof=1).over(grp),
        grp_sum=pl.col("genotype").sum().over(grp),
        grp_mean=pl.col("genotype").mean().over(grp),
        fitness_pred=((-(pl.col("phenotype") - pl.col("optimum")) ** 2 / (2 * pl.col("Vs"))).exp()
                      / (pl.col("Vs") * 2 * math.pi).sqrt()),
    )
    # copied_val. ave: mean of the others. fit / ext: the group member(s) with the best score
    # (closest to the optimum / farthest from the group mean); any tied member counts
    df = df.with_columns(
        ave_pred=(pl.col("grp_sum") - pl.col("genotype")) / (pl.col("n") - 1),
        score=pl.when(pl.col("reg") == "fit").then(-(pl.col("genotype") - pl.col("optimum")).abs())
                .otherwise((pl.col("genotype") - pl.col("grp_mean")).abs()),
    )
    df = df.with_columns(best=pl.col("score").max().over(grp))
    cands = (df.filter(close(pl.col("score"), pl.col("best"), rtol=1e-6, atol=1e-9))
             .select(*grp, pl.col("genotype").alias("cand")))
    pair = (df.select(*KEY, "group_tag", "reg", "copied_val", "ave_pred").join(cands, on=grp)
            .with_columns(hit=pl.when(pl.col("reg") == "ave")
                          .then(close(pl.col("copied_val"), pl.col("ave_pred")))
                          .otherwise(close(pl.col("copied_val"), pl.col("cand"))))
            .group_by(KEY).agg(pl.col("hit").any()))
    df = df.join(pair, on=KEY)

    checks = {
        "vars": close(pl.col("var"), pl.col("var_pred")),
        "copied_val": pl.col("hit"),
        "fitness": close(pl.col("fitness"), pl.col("fitness_pred")),
    }
    print(f"{len(ids)} sims, {df.height:,} individual-generations")
    res = df.select(pl.col("reg"), pl.col("n"), **{k: v for k, v in checks.items()})
    for k in checks:
        ok = res[k].sum()
        print(f"  {k:<11} {ok:>10,} / {res.height:,} match ({100 * ok / res.height:.3f}%)")
    by = (res.group_by("reg", "n").agg([(1 - pl.col(k).mean()).alias(f"{k}_mismatch") for k in checks])
          .sort("reg", "n"))
    bad = by.filter(pl.any_horizontal([pl.col(f"{k}_mismatch") > 0 for k in checks]))
    if bad.height:
        print("mismatch fraction by reg / n:")
        print(bad)


if __name__ == "__main__":
    main()
