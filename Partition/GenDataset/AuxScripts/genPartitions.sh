set -xe
FOLDER="Results"

jobs=5

make

#rm -rf $FOLDER
mkdir -p $FOLDER/

cd Graphs
for (( j=0; j<200; j++ ))
do
   f="GraphSamples/Graph_$j.txt"
   ./mconvert $f ConvertedGraphs/G$j.txt
   pmetis ConvertedGraphs/G$j.txt 2
   mv ConvertedGraphs/G$j.txt.part.2 Partitions/G$j.2 
done
cd ..

c=0
runCmds=()
for (( j=0; j<200; j++ ))
do
   if [ ! -f $FOLDER/G$j.out ]; then
       runCmds[c]="./prog --inputFile Graphs/ConvertedGraphs/G$j.txt --randomP 0.45 --outputFile $FOLDER/G$j.out --partitionFile Graphs/Partitions/G$j.2 --timeout 7200 &"
       c=$(($c+1))
   fi
done

i=0
while (( i<${#runCmds[@]} ))
do
  minV=$((${#runCmds[@]}-$i<$jobs ? ${#runCmds[@]}-$i : $jobs))

  pi=0
  for (( k=0; k<$minV; k++ ))
  do
    rcommand="${runCmds[i+k]}"
    eval $rcommand

    pids[${pi}]=$!
    pi=$((pi+1))
  done

  #Wait all processes
  for pid in ${pids[*]}; do
    wait $pid
  done

  i=$(($i+$jobs))
done
