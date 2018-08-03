#!/bin/bash
#
#  Compile the program with g++, run it and plot results with python
#
echo "g++ compiler and python 3 should be installed for this script to work"
compiler=`which g++`
py=`which python`
echo "compiling for sigle processor"
$compiler -o heat -O3 heat_equation.cpp 
echo "define input parameters"
$py gimmick.py

echo "running ..."
./heat
$py heatPlot.py &

#
#  Discard the executable file.
#
rm heat
#
echo "Program output written to solution.dat and Heated tile.pdf"
