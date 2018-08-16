#!/usr/bin/env bash

rm -rf build
mkdir build
cd build
cmake ../../
make

groundTruth="../../../data/MovieLens/result.txt"
QFilePath="../../../data/MovieLens/q.txt"
PFilePath="../../../data/MovieLens/p.txt"

#groundTruth="../../../data/Netflix_UM/avg_result.txt"
#QFilePath="../../../data/Netflix_UM/q-50-avg.txt"
#PFilePath="../../../data/Netflix_UM/p-50-avg.txt"

outputFilePath="result.txt"
k="10"
s="1024"


./runADiamondSampling --groundTruth $groundTruth --q_file $QFilePath --p_file $PFilePath --outputFilePath $outputFilePath --k $k --s $s