#!/bin/bash
folder=tsv_output/1t

c=0

for i in $(seq 0 50);
do
    for psi in -0.9 -0.7 -0.5 -0.3 -0.1 0 0.1 0.3 0.5 0.7 0.9;
    do
        echo "i=$i psi=$psi"
        # slim -d ID=$i -d "folder='$folder'" -d "psi11=$psi" 1t_pair.slim &> log/1t_pair_log${i}.txt &
        slim -d ID=$i -d "folder='$folder'" -d "psi11=$psi" 1t_cutoff.slim &> log/1t_pairc_log${i}.txt &
        slim -d ID=$i -d "folder='$folder'" -d "psi11=$psi" 1t_cutoff_norec.slim &> log/1t_paircnr_log${i}.txt &
        c=$((c + 2))
        if [ "$c" -eq 20 ]; then
            echo "Waiting for all jobs to finish..."
            wait
            c=0  # reset the counter if needed
        fi
    done
done