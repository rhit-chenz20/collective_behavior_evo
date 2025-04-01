

for i in $(seq 1 10);
do
    slim  -d ID=$i 1_trait.slim > log/log${i}.txt &
done



