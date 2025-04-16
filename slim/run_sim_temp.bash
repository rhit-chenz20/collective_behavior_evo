
end_gen=50

for i in $(seq 1 5);
do
    for psi in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 -0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9;
    do
        for qtl_c in 10 20 30 40 50;
        do
            echo "qtl_c: $qtl_c, psi: $psi, ID: $i"
            slim -d ID=$i -d C=$qtl_c -d psi11=$psi -d end_gen=$end_gen 1_trait.slim > log/log${i}.txt 
        done
    done
done
