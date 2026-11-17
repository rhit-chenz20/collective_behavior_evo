#!/bin/bash
# Sims that SLiM aborted because the population fell so far behind the optimum that every
# individual's fitness underflowed to 0 ("total fitness of all individuals is <= 0.0"; with
# Vs = 2000 that happens once every phenotype is more than ~1690 from the optimum).
# slurm_job.sh deletes such a sim's tsv outputs, so it is missing from the parquet data; only
# its SLiM log survives, in raw/<task>.tar.zst. This reads those logs and writes one row per
# aborted sim:
#   shift_speed reg n psi rep ID  fail_tick (first tick below 10% of the fitness at the optimum,
#   with the lag then)  crash_tick (the last 20-tick recording before the abort; the abort is
#   within the next 20 ticks)
# Run from the repo root (reads lab storage, ~1 min):
#   bash analysis/steady_tracking/server/extract_crashed_sims.sh [dataset]   (default steady_tracking)
# Output: figure/r_data/<dataset>/crashed_sims.csv
set -euo pipefail
ds=${1:-steady_tracking}
S=/fs/cbsubscb10/storage/$(id -un)/$ds
ZSTD=$HOME/miniconda3/envs/simdata/bin/zstd
out=figure/r_data/$ds/crashed_sims.csv
mkdir -p "$(dirname "$out")"
echo "shift_speed,reg,n,psi,rep,ID,fail_tick,fail_lag,crash_tick" > "$out"
for d in "$S"/summary/*/; do
    task=$(basename "$d")
    [[ -s $d/failed_sims.txt ]] || continue
    grep -q . "$d/failed_sims.txt" || continue
    arc=$S/raw/$task.tar.zst
    [[ -f $arc ]] || { echo "no archive for $task" >&2; continue; }
    tmp=$(mktemp -d)
    $ZSTD -dc "$arc" | tar -x -C "$tmp" ./log 2>/dev/null || $ZSTD -dc "$arc" | tar -x -C "$tmp"
    while read -r tag; do
        [[ -n $tag ]] || continue
        f=$tmp/log/$tag.txt
        [[ -f $f ]] || { echo "no log for $tag" >&2; continue; }
        grep -q 'total fitness of all individuals is <= 0.0' "$f" || { echo "$tag failed for another reason" >&2; continue; }
        # tag: n_<n>_psi_<psi>_sp_<sp>_reg_<reg>_<rep>_<ID>
        IFS=_ read -r _ n _ psi _ sp _ reg rep id <<< "$tag"
        fl=$(grep -m1 'failed the fitness check' "$f" | sed -E 's/^"tick ([0-9]+):.*lag ([-0-9.e]+):.*/\1,\2/' || true)
        [[ -n $fl ]] || fl=","
        ct=$(grep 'recoding at cycle' "$f" | tail -1 | grep -oE '[0-9]+' || true)
        echo "$sp,$reg,$n,$psi,$rep,$id,$fl,$ct" >> "$out"
    done < "$d/failed_sims.txt"
    rm -rf "$tmp"
done
echo "$(($(wc -l < "$out") - 1)) aborted sims written to $out"
