#!/bin/bash

gfortran -fopenmp -O3 $1.f90 -o $1.exe -L/usr/lib/x86_64-linux-gnu -lblas -L/usr/lib/x86_64-linux-gnu -llapack
sleep 2 # wait 2 seconds
 ./$1.exe &
