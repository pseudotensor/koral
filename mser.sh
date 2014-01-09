#!/bin/bash
#compiles the code using serial openmp
rm Makefile
ln -s Makefile.gcc.lap Makefile
make clean
make -j






