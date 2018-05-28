#!/usr/bin/env bash

rm -rf build
mkdir build
cd build
cmake ../../
make

groundTruth="../../../data/MovieLens/result.txt"
QFilePath="../../../data/MovieLens/q.txt"
PFilePath="../../../data/MovieLens/p.txt"
outputFilePath="result.txt"
k="10"
s="10000"


./runWedgeSampling --groundTruth $groundTruth --q_file $QFilePath --p_file $PFilePath --outputFilePath $outputFilePath --k $k  --s $s