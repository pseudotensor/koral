#!/bin/bash
#compiles the code using serial openmp
make clean
rm ko
make -j ko SERIAL=0






