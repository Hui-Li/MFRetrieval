#!/usr/bin/env bash

rm -rf build
mkdir build
cd build
cmake ../
make

#declare -a  dataset_list=("MovieLens" "Yelp" "Netflix" "KDD")
#declare -a  Q_list=("../../data/MovieLens/q.txt" "../../data/Yelp/q.txt" "../../data/Netflix/q.txt" "../../data/KDD/q.txt")
#declare -a  P_list=("../../data/MovieLens/p.txt" "../../data/Yelp/p.txt" "../../data/Netflix/p.txt" "../../data/KDD/p.txt")
#numOfDataSets=${#dataset_list[@]}

#k_list="1 2 5 10 50"
#
#SIGMA_list="0.5 0.6 0.7 0.8 0.9"
#
#scalingValue_list="1 10 100 127 1000 10000"

declare -a  dataset_list=("KDD")
declare -a  Q_list=("../../data/KDD/q.txt")
declare -a  P_list=("../../data/KDD/p.txt")
numOfDataSets=${#dataset_list[@]}

k_list="1"

SIGMA_list="0.8"

scalingValue_list="127"

################## Naive Method ####################
#for (( index=0; index < numOfDataSets; index++ )); do
#    for k in $k_list
#    do
#    drop_caches
#    ./runFEXIPRO --alg Naive --k $k --dataset ${dataset_list[$index]} --q ${Q_list[$index]} --p ${P_list[$index]}
#    done
#done
################## Naive Method ####################

################# Ball Tree ####################
#for (( index=0; index < numOfDataSets; index++ )); do
#    for k in $k_list
#    do
#    drop_caches
#    ./runFEXIPRO --alg BallTree --dataset ${dataset_list[$index]} --k $k --q ${Q_list[$index]} --p ${P_list[$index]}
#    done
#done
################# Ball Tree ####################

################# FastMKS ####################
#for (( index=0; index < numOfDataSets; index++ )); do
#    for k in $k_list
#    do
#    drop_caches
#    ./runFEXIPRO --alg FastMKS --dataset ${dataset_list[$index]} --k $k --q ${Q_list[$index]} --p ${P_list[$index]}
#    done
#done
################# FastMKS ####################

################### FEXIPRO-SIR ####################
for (( index=0; index < numOfDataSets; index++ )); do
    for k in $k_list
    do
        for SIGMA in $SIGMA_list
        do
        drop_caches
        ./runFEXIPRO --alg SIR --k $k --dataset ${dataset_list[$index]} --q ${Q_list[$index]} --p ${P_list[$index]} --SIGMA $SIGMA
        done
    done
done
################### FEXIPRO-SIR ####################

################## FEXIPRO-SR: SVD + Transformation + Incremental Prune (incr on original vector, then transformed vector) ####################
#for (( index=0; index < numOfDataSets; index++ )); do
#    for k in $k_list
#    do
#        for SIGMA in $SIGMA_list
#        do
#        drop_caches
#        ./runFEXIPRO --alg SR --k $k --dataset ${dataset_list[$index]} --q ${Q_list[$index]} --p ${P_list[$index]} --SIGMA $SIGMA
#        done
#    done
#done
################## Transformation + Incremental Prune ####################

################### FEXIPRO-S ####################
#for (( index=0; index < numOfDataSets; index++ )); do
#    for k in $k_list
#    do
#        for SIGMA in $SIGMA_list
#        do
#        drop_caches
#        ./runFEXIPRO --alg S --k $k --dataset ${dataset_list[$index]} --q ${Q_list[$index]} --p ${P_list[$index]} --SIGMA $SIGMA
#        done
#    done
#done
################## FEXIPRO-S ####################

################## FEXIPRO-I ####################
#for (( index=0; index < numOfDataSets; index++ )); do
#    for k in $k_list
#    do
#        for scalingValue in $scalingValue_list
#        do
#        drop_caches
#        ./runFEXIPRO --alg I --k $k --dataset ${dataset_list[$index]} --q ${Q_list[$index]} --p ${P_list[$index]} --scalingValue $scalingValue
#        done
#    done
#done
################## FEXIPRO-I ####################

################## FEXIPRO-SI ####################
for (( index=0; index < numOfDataSets; index++ )); do
    for k in $k_list
    do
        for SIGMA in $SIGMA_list
        do
            for scalingValue in $scalingValue_list
            do
                drop_caches
                ./runFEXIPRO --alg SI --k $k --dataset ${dataset_list[$index]} --q ${Q_list[$index]} --p ${P_list[$index]} --SIGMA $SIGMA --scalingValue $scalingValue
            done
        done
    done
done
################## FEXIPRO-SI ####################
