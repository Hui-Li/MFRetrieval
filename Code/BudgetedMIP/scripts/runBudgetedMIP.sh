#!/usr/bin/env bash

rm -rf build
mkdir build
cd build
cmake ../../
make

groundTruth="../../../data/Netflix_UM/result.txt"
QFilePath="../../../data/Netflix_UM/q-50.txt"
PFilePath="../../../data/Netflix_UM/p-50.txt"
outputFilePath="result.txt"
k="10"
budget="1024"

./runBudgetedMIP --b $budget --groundTruth $groundTruth --q_file $QFilePath --p_file $PFilePath --outputFilePath $outputFilePath --k $k