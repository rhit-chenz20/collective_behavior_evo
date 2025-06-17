#!/bin/bash
folder=tsv_output/1t

c=0

for i in $(seq 0 200);
do
    echo "i=$i"
    slim -d ID=$i -d "folder='$folder'" base.slim &> log/base_log_${i}.txt &
    slim -d ID=$i -d "folder='$folder'" base_cutoff.slim &> log/base_logc_${i}.txt &
    c=$((c + 2))
    if [ "$c" -eq 30 ]; then
        echo "Waiting for all jobs to finish..."
        wait
        c=0  # reset the counter if needed
    fi
done