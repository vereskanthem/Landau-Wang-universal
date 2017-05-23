#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "test_g/DoS-L=8/graphs/336466585.jpg"
set grid x y
set xlabel "i"
set ylabel "G(i)"
plot "test_g/DoS-L=8/336466585.dat" using 1:3 title "landau-wang-8-iteration-336466585" with lines lt rgb "red"