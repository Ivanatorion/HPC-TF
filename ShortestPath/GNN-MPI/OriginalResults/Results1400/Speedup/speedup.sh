processors=("1" "2" "4" "8" "16" "32" "48" "64" "80")

data=$(cat processTime.csv)
data=($data)

echo "0,0.00" > speedup.csv

for(( k=0; k<${#data[@]}; k++ ))
do

base=(${data[0]//,/ })
base=${base[1]}

val=(${data[k]//,/ })
val=${val[1]}

spdup=$(echo "scale=2; $base / $val" | bc)
echo "${processors[k]},$spdup" >> speedup.csv

done

gnuplot script.gnuplot
epstopdf output.eps
rm output.eps
mv output.pdf Speedup.pdf
