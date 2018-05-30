#!/usr/bin/env bash

rm -rf build-WedgePre
mkdir build-WedgePre
cd build-WedgePre
cmake ../../
make

groundTruth="../../../data/MovieLens/result.txt"
QFilePath="../../../data/MovieLens/q.txt"
PFilePath="../../../data/MovieLens/p.txt"
outputFilePath="result.txt"
k="10"
s="100"
verify="true"

./runWedgePreSampling --verify $verify --groundTruth $groundTruth --q_file $QFilePath --p_file $PFilePath --outputFilePath $outputFilePath --k $k  --s $s