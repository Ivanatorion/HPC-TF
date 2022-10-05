#!/bin/bash
#SBATCH --nodes=3
#SBATCH --ntasks=48
#SBATCH --time=23:59:00
#SBATCH --partition=draco
#SBATCH --output=../script48.out
#SBATCH --error=../script48.err

cd ..

OUTFILE="out48.txt"

MACHINEFILE="nodes.$SLURM_JOB_ID"
srun -l /bin/hostname | sort -n | awk '{print $2}' > $MACHINEFILE

make gnnmpi-nopartition &> /dev/null

graphN=100
NP=$(cat $MACHINEFILE | wc -l)
runsPerGraph=5

rm -f $OUTFILE

for (( k=0; k<$graphN; k++ ))
do
    if [ $NP != 1 ]; then
        cp Partitions/GraphPartitions/Graph_$k/Graph_$k.txt.part.$NP Partitions/Graph.part.$NP
    fi
    for (( r=0; r<$runsPerGraph; r++ ))
    do
        #echo "Graph $(($k+1)) / $graphN - Run $(($r+1)) / $runsPerGraph..."
        mpirun --mca btl ^openib --mca btl_tcp_if_include eno1 --bind-to none --mca btl_base_warn_component_unused 0 -machinefile $MACHINEFILE -np $NP gnnmpi-nopartition GNN/GraphSamples/Graph_$k.txt >> $OUTFILE
    done
    if [ $NP != 1 ]; then
        rm Partitions/Graph.part.$NP
    fi
done

total=0
count=0
lines=$(cat $OUTFILE)
for i in $lines
do
    total=$(echo $total+$i | bc )
    (( count=$count+1 ))
done

echo "scale=2; $total / $count" | bc

rm $MACHINEFILE

cd SlurmScripts
