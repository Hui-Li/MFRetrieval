#!/usr/bin/env bash

mkdir build
cd build
cmake ../
make

declare -a  dataset_list=("MovieLens" "Yelp" "Netflix" "KDD")
declare -a  Q_list=("../../data/MovieLens/q.txt" "../../data/Yelp/q.txt" "../../data/Netflix/q.txt" "../../data/KDD/q.txt")
declare -a  P_list=("../../data/MovieLens/p.txt" "../../data/Yelp/p.txt" "../../data/Netflix/p.txt" "../../data/KDD/p.txt")
numOfDataSets=${#dataset_list[@]}

k_list="1"
batchSize="1 10 100 1000 10000"

thread="0 1"

################## Naive Method ####################
#for (( index=0; index < numOfDataSets; index++ )); do
#    for k in $k_list; do
#        for s in $thread; do
#            drop_caches
#            ./main --alg MKLNaive --k $k --dataset ${dataset_list[$index]} --q ${Q_list[$index]} --p ${P_list[$index]} --singleThread $s
#        done
#    done
#done
################## Naive Method ####################

################## Batch Naive Method ####################
for (( index=0; index < numOfDataSets; index++ )); do
    for k in $k_list; do
        for b in $batchSize; do
            for s in $thread; do
                drop_caches
                ./main --alg MKLBatchNaive --batchSize $b --k $k --dataset ${dataset_list[$index]} --q ${Q_list[$index]} --p ${P_list[$index]} --singleThread $s
            done
        done
    done
done
################## Batch Naive Method ####################
