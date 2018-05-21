#!/usr/bin/env bash

rm -rf build
mkdir build
cd build
cmake ../../
make

QFilePath="../../../data/MovieLens/q.txt"
PFilePath="../../../data/MovieLens/p.txt"

outputFilePath="Netflix_UM_result.txt"
k="10"
n="-1"

./runNaive --q_file $QFilePath --p_file $PFilePath --outputFilePath $outputFilePath --k $k  --n $n