set datafile separator ','

set ylabel "Speedup"
set xlabel "Num. Processes"

set xtics 5
set ytics 1

set grid

set term eps
set output "output.eps"

plot 'speedup.csv' using 1:2 with lines notitle lw 5 lt rgb "red"

