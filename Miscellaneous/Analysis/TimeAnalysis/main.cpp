#include "util/Base.h"
#include "structs/Matrix.h"
#include "alg/BasicNormPrune.h"
#include "alg/BasicNormPruneIncr.h"

void scanAndCompute(string qDataPath, string pDataPath, int k) {
    Matrix q;
    Matrix p;
    q.readData(qDataPath);
    p.readData(pDataPath);

    cout << "q:" << q.rowNum << endl;
    cout << "p:" << p.rowNum << endl;

    // step 1: record when to stop
    unsigned long long computationAmount;
    unsigned long long heapOperationAmount1;
    unsigned long long heapOperationAmount2;
    basicNormPruneRecord(k, q, p, computationAmount, heapOperationAmount1, heapOperationAmount2);

    cout << "computationAmount: " << computationAmount << endl;
    cout << "heapOperationAmount1: " << heapOperationAmount1 << endl;
    cout << "heapOperationAmount2: " << heapOperationAmount2 << endl;

    cout << "step 1 finish" << endl;

    // step 2: only scan
    basicNormPruneOnlyScan(k, q, p, computationAmount, heapOperationAmount1, heapOperationAmount2);

//    basicNormPruneOnlyCompute(q.rowNum, q.colNum, readFlag);

    cout << "step 2 finish" << endl;

}

void orderTimeForIncr(string qDataPath, string pDataPath, int k) {

    Matrix q;
    Matrix p;

    q.readData(qDataPath);
    p.readData(pDataPath);

    cout << "q:" << q.rowNum << endl;
    cout << "p:" << p.rowNum << endl;

    basicNormPruneIncr(k, 10, q, p);

}


int main(int argc, char **argv) {

    string qDataPath = "./q.txt";
    string pDataPath = "./p.txt";

    int k = 1;

    scanAndCompute(qDataPath, pDataPath, k);

    // not useable. Order by q requires recomupting subnorm!!!
    orderTimeForIncr(qDataPath, pDataPath, k);

    return 0;
}
