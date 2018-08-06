#!/bin/bash
#
#  Compile the program with G++
#
echo "g++ compiler should be installed for this script to work"
g++ -fopenmp heat_equation.cpp -lm -o heat_p
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./heat_p > output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./heat_p >> output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./heat_p >> output.txt
#
#  Discard the executable file.
#
rm heat_p
#
echo "Program output written to output.txt"
