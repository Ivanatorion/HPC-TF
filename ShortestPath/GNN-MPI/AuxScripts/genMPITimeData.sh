set -xe
cores=("1" "2" "3" "4" "5" "6")
export RST_BUFFER_SIZE=80000000
FOLDER="Results/PartitionComm"
NUM_GRAPHS=50
cPath=`realpath .`

make gnnmpi
make AuxPrograms/mpitimecalc
cp GNN/model_weights_low.dat GNN/model_weights.dat
cp -r GNN/GraphSamples1400-1500 GNN/GraphSamples

rm -rf $FOLDER

for (( k=0; k<${#cores[@]}; k++ ))
do
    mkdir -p $FOLDER/Metis/Cores${cores[k]}
    mkdir -p $FOLDER/SimplePartition/Cores${cores[k]}
done

echo Graph,Processes,CutEdges > $FOLDER/Metis/cutEdges.csv
echo Graph,Processes,CutEdges > $FOLDER/SimplePartition/cutEdges.csv

for (( j=0; j<$NUM_GRAPHS; j++ ))
do
   f="GNN/GraphSamples/Graph_$j.txt"
   for (( k=0; k<${#cores[@]}; k++ ))
   do
       mpirun -np ${cores[k]} gnnmpi $f 1 > outtimedata.txt
       AuxScripts/processRastro.sh
       data=$(cat outtimedata.txt)
       data=($data)
       base=${data[5]}
       
       echo $j,${cores[k]},${data[5]} >> $FOLDER/Metis/cutEdges.csv
       rm outtimedata.txt
       mv mpitime.csv $FOLDER/Metis/Cores${cores[k]}/G$j.csv
    done
done

for (( j=0; j<$NUM_GRAPHS; j++ ))
do
    f="GNN/GraphSamples/Graph_$j.txt"
    for (( k=0; k<${#cores[@]}; k++ ))
    do
        mpirun -np ${cores[k]} gnnmpi $f 0 > outtimedata.txt
        AuxScripts/processRastro.sh
        data=$(cat outtimedata.txt)
        data=($data)
        base=${data[5]}
        
        echo $j,${cores[k]},${data[5]} >> $FOLDER/SimplePartition/cutEdges.csv
        rm outtimedata.txt
        mv mpitime.csv $FOLDER/SimplePartition/Cores${cores[k]}/G$j.csv
     done
done

rm GNN/model_weights.dat
rm -r GNN/GraphSamples

echo Cores,METISCutEdges,METISTime,METISCommRatio,SIMPLECutEdges,SIMPLETime,SIMPLECommRatio > "$FOLDER/avg.csv"

set +x
for (( k=0; k<${#cores[@]}; k++ ))
do
    METIScutEdgeAvg=0
    f="$FOLDER/Metis/cutEdges.csv"
    edgedata=$(cat $f)
    edgedata=$(awk -F',' '{ for( i=1; i<=NF; i++ ) print $i }' <<< "$edgedata")
    edgedata=($edgedata)

    METISsumCommRatio=0.0000
    METISsumTime=0.0000
    for (( j=0; j<$NUM_GRAPHS; j++ ))
    do
        f="$FOLDER/Metis/Cores${cores[k]}/G$j.csv"
        data=$(cat $f)
        data=$(awk -F',' '{ for( i=1; i<=NF; i++ ) print $i }' <<< "$data")
        data=($data)
        cRatio=0.00
        for (( rankIdx=0; rankIdx<${cores[k]}; rankIdx++ ))
        do
            cRatio=$(awk "BEGIN{ print $cRatio + ${data[$((9 + rankIdx * 5))]} }")
        done
        cRatio=$(awk "BEGIN{ print $cRatio / ${cores[k]} }")
        METISsumCommRatio=$(awk "BEGIN{ print $METISsumCommRatio + $cRatio }")

        METISsumTime=$(awk "BEGIN{ print $METISsumTime + ${data[7]} }")
        METIScutEdgeAvg=$(awk "BEGIN{ print $METIScutEdgeAvg + ${edgedata[$((5 + 3 * k + j * 18))]} }")
    done
    METISsumCommRatio=$(awk "BEGIN{ print $METISsumCommRatio / $NUM_GRAPHS.0 }")
    METISsumTime=$(awk "BEGIN{ print $METISsumTime / $NUM_GRAPHS.0 }")
    METIScutEdgeAvg=$(awk "BEGIN{ print $METIScutEdgeAvg / $NUM_GRAPHS.0 }")

    SIMPLEcutEdgeAvg=0
    f="$FOLDER/SimplePartition/cutEdges.csv"
    edgedata=$(cat $f)
    edgedata=$(awk -F',' '{ for( i=1; i<=NF; i++ ) print $i }' <<< "$edgedata")
    edgedata=($edgedata)

    SIMPLEsumCommRatio=0.0000
    SIMPLEsumTime=0.0000
    for (( j=0; j<$NUM_GRAPHS; j++ ))
    do
        f="$FOLDER/SimplePartition/Cores${cores[k]}/G$j.csv"
        data=$(cat $f)
        data=$(awk -F',' '{ for( i=1; i<=NF; i++ ) print $i }' <<< "$data")
        data=($data)
        cRatio=0.00
        for (( rankIdx=0; rankIdx<${cores[k]}; rankIdx++ ))
        do
            cRatio=$(awk "BEGIN{ print $cRatio + ${data[$((9 + rankIdx * 5))]} }")
        done
        cRatio=$(awk "BEGIN{ print $cRatio / ${cores[k]} }")
        SIMPLEsumCommRatio=$(awk "BEGIN{ print $SIMPLEsumCommRatio + $cRatio }")

        SIMPLEsumTime=$(awk "BEGIN{ print $SIMPLEsumTime + ${data[7]} }")
        SIMPLEcutEdgeAvg=$(awk "BEGIN{ print $SIMPLEcutEdgeAvg + ${edgedata[$((5 + 3 * k + j * 18))]} }")
    done
    SIMPLEsumCommRatio=$(awk "BEGIN{ print $SIMPLEsumCommRatio / $NUM_GRAPHS.0 }")
    SIMPLEsumTime=$(awk "BEGIN{ print $SIMPLEsumTime / $NUM_GRAPHS.0 }")
    SIMPLEcutEdgeAvg=$(awk "BEGIN{ print $SIMPLEcutEdgeAvg / $NUM_GRAPHS.0 }")

    echo ${cores[k]},$METIScutEdgeAvg,$METISsumTime,$METISsumCommRatio,$SIMPLEcutEdgeAvg,$SIMPLEsumTime,$SIMPLEsumCommRatio >> "$FOLDER/avg.csv"
done
set -x

cp AuxScripts/partition.gnuplot $FOLDER/script.gnuplot
cd $FOLDER

gnuplot script.gnuplot
epstopdf Partition.eps
rm Partition.eps

cd $cPath
