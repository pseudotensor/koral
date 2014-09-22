#!/bin/bash
#compiles the code using serial openmp
make clean
rm ana avg ko
make -j SERIAL=1






