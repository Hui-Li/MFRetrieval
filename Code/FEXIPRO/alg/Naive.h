#ifndef NAIVE_H
#define NAIVE_H

#include "../util/Base.h"
#include "../util/FileUtil.h"
#include "../structs/Matrix.h"
#include "../util/Calculator.h"
#include "../util/Monitor.h"

/**
 * Naive solution.
 */
void naive(const int k, const Matrix &q, const Matrix &p) {

    Monitor t;
    t.start();

    vector<vector<VectorElement> > results(q.rowNum, vector<VectorElement>());

    for (int qID = 0; qID < q.rowNum; qID++) {

        const double *qRow = q.getRowPtr(qID);

        vector<VectorElement> &heap = results[qID];
        heap.resize(k);

        for (int pID = 0; pID < k; pID++) {
            const double *pRow = p.getRowPtr(pID);
            double innerProduct = Calculator::innerProduct(qRow, pRow, p.colNum);
            heap[pID] = VectorElement(pID, innerProduct);
        }

        make_heap(heap.begin(), heap.end(), greater<VectorElement>());
        double lowerBound = heap.front().data;

        for (int pID = k; pID < p.rowNum; pID++) {
            const double *pRow = p.getRowPtr(pID);
            double innerProduct = Calculator::innerProduct(qRow, pRow, p.colNum);

            if (innerProduct > lowerBound) {
                pop_heap(heap.begin(), heap.end(), greater<VectorElement>());
                heap.pop_back();
                heap.push_back(VectorElement(pID, innerProduct));
                push_heap(heap.begin(), heap.end(), greater<VectorElement>());
                lowerBound = heap.front().data;
            }
        }
    }

    t.stop();

    Logger::Log("time: " + to_string(t.getElapsedTime()) + " secs");

    if (Conf::outputResult) {
        string resultFileName = Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(k) + ".txt";
        FileUtil::outputResult(k, results, resultFileName);
    }
}


#endif //NAIVE_H
