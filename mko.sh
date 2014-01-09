#!/bin/bash
#compiles the code using mpi
rm Makefile
ln -s Makefile.mpigcc.lap Makefile
make clean
make ko -j






