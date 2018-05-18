#!/usr/bin/env bash

rm -rf build
mkdir build
cd build
cmake ../../
make

compareWithNaive="true"
QFilePath="../../../data/MovieLens/q.txt"
PFilePath="../../../data/MovieLens/p.txt"
outputFilePath="result.txt"
k="10"
s="1000"


./runWedgeSampling --compareWithNaive $compareWithNaive --q_file $QFilePath --p_file $PFilePath --outputFilePath $outputFilePath --k $k  --s $s