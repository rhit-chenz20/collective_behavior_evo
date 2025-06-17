#!/bin/bash
folder=tsv_output/newM

c=0

for i in $(seq 0 9);
do
    
    echo "i=$i"
    # slim -d ID=$i -d "folder='$folder'" -d "psi11=$psi" 1t_pair.slim &> log/1t_pair_log${i}.txt &
    slim -d ID=$i -d "folder='$folder'" newmodel.slim &> log/newM_${i}.txt &
    c=$((c + 1))
    if [ "$c" -eq 10 ]; then
        echo "Waiting for all jobs to finish..."
        wait
        c=0  # reset the counter if needed
    fi

done