#!/bin/bash
# Run this from the cluster login node (cbsulogin), inside this directory under /home.
# Resubmits every task of the full run (submit_slurm.sh) that has no $OUT/summary/<chunk>/_DONE:
# timed out, failed, or never ran. A rerun task first removes its chunk's old outputs
# (slurm_job.sh step 4), so nothing needs cleaning by hand.
#
# Reruns use the same EXCLUDE_NODES and TIME_LIMIT as submit_slurm.sh (pipeline_config.sh).
#
#   bash rerun_failed.sh            list the tasks to rerun, submit nothing
#   bash rerun_failed.sh --submit   submit them

set -euo pipefail

MAX_SPEC=${MAX_SPEC:-1000}   # longest --array list per sbatch call

source pipeline_config.sh
mount_lab_storage > /dev/null

# same task -> chunk name mapping as slurm_job.sh
n_gap=${#GAPS[@]}; n_sz=${#SZS[@]}; n_n=${#NS[@]}; n_psi=${#PSIS[@]}; n_reg=${#REGS[@]}
todo=()
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
    [[ -f $OUT/summary/${part}/_DONE ]] || todo+=($task)
done
echo "${#todo[@]} of ${N_TASKS} tasks have no _DONE in ${OUT}/summary"
(( ${#todo[@]} )) || exit 0

# array indices must stay below MAX_ARRAY, so keep submit_slurm.sh's offsets and give each
# array the list of its tasks that need a rerun
submit=0; [[ ${1:-} == "--submit" ]] && submit=1
(( submit )) && mkdir -p /home/$USER/slurm_logs
for (( off=0; off<N_TASKS; off+=MAX_ARRAY )); do
    ids=()
    for t in "${todo[@]}"; do
        (( t >= off && t < off + MAX_ARRAY )) && ids+=($(( t - off )))
    done
    (( ${#ids[@]} )) || continue
    # collapse runs of consecutive ids into ranges (370,371,372 -> 370-372), then split them
    # into --array lists of at most MAX_SPEC characters: sbatch rejects long ones ("Pathname
    # of a file, directory or other parameter too long"), so one offset may need several calls
    mapfile -t lists < <(printf '%s\n' "${ids[@]}" | awk -v max=$MAX_SPEC '
        function add(r) { if (cur != "" && length(cur) + 1 + length(r) > max) { print cur; cur = "" }
                          cur = (cur == "" ? r : cur "," r) }
        NR==1{s=p=$1; next}
        $1==p+1{p=$1; next}
        {add(p>s ? s "-" p : s); s=p=$1}
        END{add(p>s ? s "-" p : s); print cur}')
    for list in "${lists[@]}"; do
        n_list=$(echo "$list" | tr ',' '\n' | awk -F- '{n += (NF>1 ? $2-$1+1 : 1)} END{print n}')
        if (( submit )); then
            id=$(sbatch --parsable --time=$TIME_LIMIT --exclude="$EXCLUDE_NODES" --array=$list \
                --export=ALL,TASK_OFFSET=${off} slurm_job.sh)
            echo "sims:    ${id}  (${n_list} tasks, offset ${off})"
        else
            echo "offset ${off}: ${n_list} tasks, --array of ${#list} characters"
        fi
    done
done
if (( submit )); then
    echo "finished tasks: ls ${OUT}/summary/*/_DONE | wc -l   (expect ${N_TASKS})"
else
    echo "nothing submitted; to submit (--time=$TIME_LIMIT --exclude=$EXCLUDE_NODES): bash rerun_failed.sh --submit"
fi
