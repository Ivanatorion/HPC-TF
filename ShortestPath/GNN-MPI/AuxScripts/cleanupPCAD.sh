process=("1" "2" "4" "8" "16" "32" "48" "64" "80")

rm -f *.err
rm -f processTime.csv

for (( k=0; k<${#process[@]}; k++ ))
do
  p=${process[k]}
  times=$(cat script$p.out)
  echo "$p,$times" >> processTime.csv
done
