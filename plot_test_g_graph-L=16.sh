#!/bin/bash
gnuplot test_g/DoS-L=16/temp/*.plot
convert -delay 100 -loop 0 test_g/DoS-L=16/graphs/{1..20}.jpg animate-DoS-L=16.gif
