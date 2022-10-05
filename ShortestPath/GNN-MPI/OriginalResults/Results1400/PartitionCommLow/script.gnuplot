set datafile separator ','

SIMPLE = "#1f2763"; METIS = "#308016"

set auto y

set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
set xtic scale 0

set term eps color size 6, 8
set output "Partition.eps"

set key center top outside

set multiplot layout 3,1
set ylabel "Cut Edges"
set xlabel ""
plot 'avg.csv' using 2:xtic(1) title "Metis Partition" lc rgb METIS, '' u 5 title "Simple Partition" lc rgb SIMPLE

set ylabel "Process Time"
set xlabel ""
set yrange[0:20]
plot 'avg.csv' using 3:xtic(1) notitle lc rgb METIS, '' u 6 notitle lc rgb SIMPLE

set ylabel "Comm. Time (%)"
set xlabel "Num. Processes"
set yrange[0:100]
plot 'avg.csv' using 4:xtic(1) notitle lc rgb METIS, '' u 7 notitle lc rgb SIMPLE
unset multiplot
