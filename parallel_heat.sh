#!/bin/bash
#
#  Compile the program with g++
#
echo "g++ compiler should be installed for this script to work"
g++ -fopenmp heat_equation_omp.cpp -lm -o heat_p
file="input"
if [ ! -f $file ]; then
    python gimmick.py
fi
#
#  Run with 1, 2, and 4 threads.
#
rm -f output.txt
for i in {1,2,4};
do
    echo "Run with $i threads"
    export OMP_NUM_THREADS=$i
    ./heat_p >> output.txt
done
#
#  Discard the executable file.
#
rm heat_p
#
echo "Program output written to output.txt"
