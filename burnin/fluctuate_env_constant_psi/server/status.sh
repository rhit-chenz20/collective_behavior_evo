#!/bin/bash
# Progress of the run set up in pipeline_config.sh (read-only: submits and changes nothing).
#   bash status.sh              counts of tasks and sims
#   bash status.sh --list       also list the failed tasks and failed sims
# A task (one combination x REPS_PER_TASK reps) writes $OUT/summary/<chunk>/ only at its end
# (slurm_job.sh step 4), and _DONE in it only if every step succeeded. So:
#   done     _DONE present
#   failed   summary folder but no _DONE (a sim, conversion, archive or copy failed)
#   running / pending   in squeue now
#   not started         none of the above (includes tasks killed before step 4, e.g. by
#                       the time limit: see the slurm logs)
# Rerun failed and unfinished tasks with rerun_failed.sh.

set -uo pipefail
cd "$(dirname "$0")"
source pipeline_config.sh
mount_lab_storage > /dev/null 2>&1

list=0; [[ ${1:-} == "--list" ]] && list=1

# same task -> chunk name mapping as slurm_job.sh
n_gap=${#GAPS[@]}; n_sz=${#SZS[@]}; n_n=${#NS[@]}; n_psi=${#PSIS[@]}; n_reg=${#REGS[@]}
n_done=0; n_fail=0; n_none=0; sims_done=0; sims_failed=0
failed_tasks=(); failed_sim_list=()
for (( task=0; task<N_TASKS; task++ )); do
    chunk=$(( task % N_CHUNKS ))
    r=$(( task / N_CHUNKS ))
    reg_i=$(( r % n_reg ));  r=$(( r / n_reg ))
    psi_i=$(( r % n_psi ));  r=$(( r / n_psi ))
    n_i=$(( r % n_n ));      r=$(( r / n_n ))
    sz_i=$(( r % n_sz ));    r=$(( r / n_sz ))
    gap_i=$r
    rep_first=$(( chunk * REPS_PER_TASK + 1 ))
    rep_last=$(( rep_first + REPS_PER_TASK - 1 ))
    (( rep_last > REPS )) && rep_last=$REPS
    part="gap_${GAPS[$gap_i]}_sz_${SZS[$sz_i]}_reg_${REGS[$reg_i]}_n_${NS[$n_i]}_psi_${PSIS[$psi_i]}_reps_${rep_first}-${rep_last}"
    dir=$OUT/summary/$part
    n_reps=$(( rep_last - rep_first + 1 ))
    if [[ -f $dir/_DONE ]]; then
        (( n_done++ )); (( sims_done += n_reps ))
    elif [[ -d $dir ]]; then
        (( n_fail++ )); failed_tasks+=("$task $part")
        if [[ -s $dir/failed_sims.txt ]]; then
            k=$(grep -c . "$dir/failed_sims.txt")
            (( sims_failed += k )); (( sims_done += n_reps - k ))
            while read -r s; do [[ -n $s ]] && failed_sim_list+=("$part: $s"); done < "$dir/failed_sims.txt"
        fi
    else
        (( n_none++ ))
    fi
done

# jobs of this pipeline in the queue (job name from slurm_job.sh's #SBATCH --job-name)
job_name=$(grep -m1 -oP '(?<=--job-name=)\S+' slurm_job.sh)
n_run=0; n_pend=0
if command -v squeue > /dev/null; then
    n_run=$(squeue -h -u "$USER" -t RUNNING ${job_name:+-n "$job_name"} -r 2>/dev/null | wc -l)
    n_pend=$(squeue -h -u "$USER" -t PENDING ${job_name:+-n "$job_name"} -r 2>/dev/null | wc -l)
fi
n_left=$(( N_TASKS - n_done ))
total_sims=$(( N_COMBOS * REPS ))

pct() { awk -v a=$1 -v b=$2 'BEGIN{printf "%.1f", b ? 100*a/b : 0}'; }
echo "Run: $OUT"
echo "Tasks: $N_TASKS ($N_COMBOS combinations x $N_CHUNKS chunks of $REPS_PER_TASK reps)"
printf "  done          %6d  (%s%%)\n" $n_done "$(pct $n_done $N_TASKS)"
printf "  failed        %6d  (summary but no _DONE)\n" $n_fail
printf "  running       %6d  (squeue%s)\n" $n_run "${job_name:+, job name $job_name}"
printf "  pending       %6d\n" $n_pend
printf "  not started   %6d  (no summary folder; includes running, pending and killed tasks)\n" $n_none
printf "  left to do    %6d  (all but done)\n" $n_left
echo "Sims: $total_sims"
printf "  finished      %6d  (%s%%; in done tasks, plus the successful sims of failed tasks)\n" $sims_done "$(pct $sims_done $total_sims)"
printf "  failed        %6d  (listed in failed_sims.txt)\n" $sims_failed
if (( list )); then
    (( n_fail )) && { echo; echo "Failed tasks (task index, chunk):"; printf '  %s\n' "${failed_tasks[@]}"; }
    (( sims_failed )) && { echo; echo "Failed sims:"; printf '  %s\n' "${failed_sim_list[@]}"; }
fi
(( n_fail || n_none )) && echo "Rerun failed / unfinished tasks once nothing is running: bash rerun_failed.sh (--submit)"
exit 0
