#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "graph/DensityStat_L=8/gi.jpg"
set grid x y
set xlabel "i"
set ylabel "g(i)"
plot "results/DensityStat_L=8.dat" using 1:3 title "landau-wang-8" with lines lt rgb "red"