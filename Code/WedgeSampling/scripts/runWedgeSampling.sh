#!/usr/bin/env bash

rm -rf build-Wedge
mkdir build-Wedge
cd build-Wedge
cmake ../../
make

groundTruth="../../../data/MovieLens/result.txt"
QFilePath="../../../data/MovieLens/q.txt"
PFilePath="../../../data/MovieLens/p.txt"
outputFilePath="result.txt"
k="10"
s="1024"
verify="true"

./runWedgeSampling --verify $verify --groundTruth $groundTruth --q_file $QFilePath --p_file $PFilePath --outputFilePath $outputFilePath --k $k  --s $s