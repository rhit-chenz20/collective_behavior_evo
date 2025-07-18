#!/bin/bash
folder=tsv_output/1t/social

c=0

for i in $(seq 0 50);
do
    for psi in -0.9 -0.7 -0.5 -0.3 -0.1 0 0.1 0.3 0.5 0.7 0.9;
    do
        for j in $(seq 1 2 21)
        do
            echo "i=$i psi=$psi c=$j"
            slim -d ID=$i -d "folder='$folder'" -d "K=$j" -d "psi11=$psi" 1t_cutoff_max.slim &> log/1t_max_${j}_${i}.txt &
            c=$((c + 1))
            if [ "$c" -eq 30 ]; then
                echo "Waiting for all jobs to finish..."
                wait
                c=0  # reset the counter if needed
            fi
        done
    done
done