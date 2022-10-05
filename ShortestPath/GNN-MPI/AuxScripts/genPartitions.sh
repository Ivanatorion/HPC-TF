#Generate metis partitions for use in gnnmpi-nopartition
#gnnmpi-nopartition reads the partition from a file rather than creating a partition in runtime
#Useful because it is not necessary to link against metis c library

set -xe

parts=("2" "4" "8" "16" "32" "48" "64" "80")

make

rm -rf Partitions/GraphPartitions
mkdir -p Partitions/GraphPartitions

rm -rf Partitions/ConvertedGraphs
mkdir -p Partitions/ConvertedGraphs

cd GNN/GraphSamples
for f in *.txt
do
    ../../AuxPrograms/mconvert $f ../../Partitions/ConvertedGraphs/$f 
done
cd ../../Partitions/ConvertedGraphs

for f in *.txt
do
    fileName=${f%.txt}
    mkdir ../GraphPartitions/$fileName
    for (( k=0; k<${#parts[@]}; k++ ))
    do
        nParts=${parts[k]}
        pmetis $f $nParts
    done
    mv $f.part* ../GraphPartitions/$fileName/.
done

cd ../..
