# Settings for run_sim.sh and run_combo.sh on the Mac (the cluster has its own in ../server).
# Values already set in the environment (e.g. DIR_NAME=test_run via --export) take precedence.

DIR_NAME=${DIR_NAME:-cons_env_moving_opt}

# burnin pop files (input): cons_n_<n>_psi_<psi>_reg_<reg>_<rep>.pop, one per n/psi/reg/rep;
# every shift speed starts from the same burnin
# Everything stays under <repo>/data (git-ignored) on the local disk; set DATA_ROOT to use an
# external drive (no spaces in the path).
DATA_ROOT=${DATA_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)/data}
OUT=${OUT:-${DATA_ROOT}/${DIR_NAME}}
POP_DIR=${POP_DIR:-${DATA_ROOT}/cons_env_pop}
WORK_ROOT=${OUT}/work
# parallel combos, default one per physical core
N_JOBS=${N_JOBS:-$(sysctl -n hw.physicalcpu)}
# stops macOS tar adding ._ files, which would fail the archive file count
export COPYFILE_DISABLE=1

# parameter grid: one array task per (shift speed, n, psi, reg) combination, running all
# replicates in a row. Task index order: speed > n > psi > reg (reg varies fastest).
# The optimum moves by shift_speed every tick from tick 10000 (stops at 11000 or when
# mean fitness drops below half the fitness at the optimum, see base_script.slim).
SPEEDS=(0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1)
NS=(2 4 5 8 10 20 25 40 50 100)
PSIS=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
REGS=(ave ext fit)
REPS=${REPS:-10}        # replicates per task; a rep keeps its sim ID whatever REPS is
N_COMBOS=$(( ${#SPEEDS[@]} * ${#NS[@]} * ${#PSIS[@]} * ${#REGS[@]} ))   # 3300

# tsv output folders converted to parquet; log/ holds .txt and is only archived
TSV_TYPES=(phenotype genotype data vars fitness allele copied_val group_tag metadata)
ARCHIVE_TYPES=("${TSV_TYPES[@]}" log)

PY=${PY:-python3}
ZSTD=${ZSTD:-zstd}
SLIM=${SLIM:-slim}
NTFY_TOPIC=mov_opt_sim   # run_sim.sh posts updates to ntfy.sh/$NTFY_TOPIC
