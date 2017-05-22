#!/bin/bash
gnuplot test_g/DoS-L=0/temp/*.plot
convert -delay 100 -loop 0 test_g/DoS-L=0/graphs/{1..20}.jpg animate-DoS-L=0.gif
