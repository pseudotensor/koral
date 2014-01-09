#!/bin/bash
#compiles the code using serial openmp
rm Makefile
ln -s Makefile.mpigcc.lap Makefile
make clean
make -j






