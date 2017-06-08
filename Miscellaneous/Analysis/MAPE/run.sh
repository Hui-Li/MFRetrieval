#!/usr/bin/env bash
# for server
export LD_LIBRARY_PATH=/home/hli2/libs/gmp/lib:/home/hli2/libs/mpfr/lib:/home/hli2/libs/mpc/lib

mkdir build
cd build
cmake -DCMAKE_C_COMPILER='/home/hli2/libs/gcc/bin/gcc' -DCMAKE_CXX_COMPILER='/home/hli2/libs/gcc/bin/g++' -DCMAKE_EXE_LINKER_FLAGS='-Wl,-rpath,/home/hli2/libs/gcc/lib64' ../
make

./main