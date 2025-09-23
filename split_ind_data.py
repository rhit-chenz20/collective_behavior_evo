import os
import re
from typing import Dict, IO

GEN_SUFFIX_RE = re.compile(r"_gen\d+\.tsv$", re.IGNORECASE)
SPLIT_NAME_RE = re.compile(r"^(?P<base>.+)_gen\d+\.tsv$", re.IGNORECASE)

def split_all_ind_files_by_generation(root_folder: str) -> None:
    """
    For each directory under root_folder:
      - Identify raw *.tsv files that do NOT end with _gen<INT>.tsv AND
        also do NOT have any sibling file named <base>_gen<INT>.tsv.
      - For each such raw file, split lines by the FIRST column (generation, must be int)
        and write <base>_gen<GEN>.tsv files next to it.
    """
    for dirpath, _, filenames in os.walk(root_folder):
        # keep only .tsv
        tsvs = [f for f in filenames if f.endswith(".tsv")]

        # Build a set of base names that already have split outputs in this directory
        bases_with_split = set()
        for fname in tsvs:
            m = SPLIT_NAME_RE.match(fname)
            if m:
                bases_with_split.add(m.group("base"))

        # Decide which files are raw & unsplit (eligible for processing)
        raw_candidates = []
        for fname in tsvs:
            if GEN_SUFFIX_RE.search(fname):
                # already looks like a split file -> skip
                continue
            base = os.path.splitext(fname)[0]
            if base in bases_with_split:
                # there's already <base>_genXXXXX.tsv present -> original is considered "already split"
                continue
            raw_candidates.append(fname)

        # Process each raw file
        for fname in raw_candidates:
            file_path = os.path.join(dirpath, fname)
            base = os.path.splitext(fname)[0]
            output_dir = dirpath

            open_files: Dict[int, IO] = {}  # gen -> file handle
            counts: Dict[int, int] = {}

            def get_outfile(gen_int: int) -> IO:
                f = open_files.get(gen_int)
                if f is None:
                    out_name = f"{base}_gen{gen_int}.tsv"
                    out_path = os.path.join(output_dir, out_name)
                    f = open(out_path, "w", newline="")
                    open_files[gen_int] = f
                    counts[gen_int] = 0
                return f

            total_lines = 0
            with open(file_path, "r") as infile:
                for line in infile:
                    total_lines += 1
                    fields = line.rstrip("\n").split("\t")
                    if not fields:
                        continue
                    # FIRST column is generation; skip non-integer rows (headers/malformed)
                    try:
                        gen = int(fields[0])
                    except (ValueError, TypeError):
                        continue
                    out_f = get_outfile(gen)
                    out_f.write(line)
                    counts[gen] += 1

            # Close any opened split files
            for f in open_files.values():
                f.close()

            if counts:
                gens_info = ", ".join(f"{g}({c})" for g, c in sorted(counts.items()))
                print(f"Processed {file_path}: total={total_lines}; wrote -> {gens_info}")
            else:
                print(f"Processed {file_path}: total={total_lines}; no valid generation lines found.")

# Example:
split_all_ind_files_by_generation("data/n_group_extreme/ind")
split_all_ind_files_by_generation("data/n_group_var_model/ind")
