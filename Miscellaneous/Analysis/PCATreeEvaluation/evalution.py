#-*-coding:utf-8-*-
from __future__ import division

from math import sqrt

def innerProduct(qVec, pVec):

    value = 0
    for i in range(0, len(qVec)):
        value = value + qVec[i] * pVec[i]

    return value

def calPrecision(k, ground_truth_file, pcatree_file):

    ground_truth = {}
    pcatree_result = {}

    with open(ground_truth_file) as input:
        for line in input.readlines():
            line = line.strip()
            paras = line.split(" ")
            recommend_list = ground_truth.get(paras[0], [])
            recommend_list.append(paras[1])
            ground_truth[paras[0]] = recommend_list

    with open(pcatree_file) as input:
        for line in input.readlines():
            line = line.strip()
            paras = line.split(" ")
            recommend_list = pcatree_result.get(paras[0], [])
            recommend_list.append(paras[1])
            pcatree_result[paras[0]] = recommend_list

    total = 0
    for user_id, recommend_list in ground_truth.iteritems():
        total = total + len(set(recommend_list).intersection(set(pcatree_result[user_id])))

    total = total / k / len(ground_truth)

    print ground_truth_file
    print pcatree_file
    print str(total)
    print "--------------------------------------------------------------"

def calRMSE(k, ground_truth_file, pcatree_file, q, p):

    ground_truth = {}
    pcatree_result = {}

    with open(ground_truth_file) as input:
        for line in input.readlines():
            line = line.strip()
            paras = line.split(" ")
            recommend_list = ground_truth.get(paras[0], [])
            recommend_list.append(paras[2])
            ground_truth[paras[0]] = recommend_list

    with open(pcatree_file) as input:
        for line in input.readlines():
            line = line.strip()
            paras = line.split(" ")
            recommend_list = pcatree_result.get(paras[0], [])
            qid = int(paras[0])
            pid = int(paras[1])

            recommend_list.append(innerProduct(q[qid], p[pid]))
            pcatree_result[paras[0]] = recommend_list

    rmse = 0
    for user_id, recommend_list in ground_truth.iteritems():
        recommend_list.sort()
        pcatree_result[user_id].sort()

        for i in range(k):
            tmp = float(recommend_list[i]) - float(pcatree_result[user_id][i])
            # print str(recommend_list[i]) + "," + str(pcatree_result[user_id][i])
            rmse = rmse + tmp * tmp

    rmse = sqrt(rmse / k / len(ground_truth))

    print ground_truth_file
    print pcatree_file
    print str(rmse)
    print "--------------------------------------------------------------"


if __name__ == "__main__":

    # datasets = ["MovieLens", "Yelp", "Netflix", "KDD"]
    datasets = ["Netflix", "KDD"]
    # k_list = [1, 2, 5, 10, 50]
    k_list = [10, 50]

    for dataset in datasets:

        qFile = dataset + "/q.txt"
        pFile = dataset + "/p.txt"

        q = []
        p = []

        with open(qFile) as input:
            for line in input.readlines():
                line = line.strip()
                if line == "":
                    continue

                q.append(map(float, line.split(",")))

        with open(pFile) as input:
            for line in input.readlines():
                line = line.strip()
                if line == "":
                    continue

                p.append(map(float, line.split(",")))

        for k in k_list:
            print dataset + "," + str(k)
            ground_truth_file = dataset + "/result-" + dataset + "-Naive-" + str(k) + ".txt"
            pcatree_file = dataset + "/result-" + dataset + "-PCATree-6-" + str(k) + ".txt"

            calRMSE(k, ground_truth_file, pcatree_file, q, p)
            # calPrecision(k, ground_truth_file, pcatree_file)