#!/usr/bin/env bash

rm -rf build
mkdir build
cd build
cmake ../../
make

groundTruth="../../../data/MovieLens/result.txt"
QFilePath="../../../data/MovieLens/q.txt"
PFilePath="../../../data/MovieLens/p.txt"

verify="true"
outputFilePath="result.txt"
k="10"
search_k="512"
s="512"


./runDiamondSampling --verify $verify --groundTruth $groundTruth --q_file $QFilePath --p_file $PFilePath --outputFilePath $outputFilePath --k $k  --search_k $search_k --s $s