#ifndef BASICNORMPRUNE_H
#define BASICNORMPRUNE_H

#include "../util/Base.h"
#include "../util/FileUtil.h"
#include "../util/Calculator.h"
#include "../util/Monitor.h"
#include "../structs/VectorElement.h"
#include "../structs/Matrix.h"
#include "../structs/ExtendMatrix.h"
#include "../structs/ExtendMatrixRow.h"


void basicNormPruneRecord(const int k, Matrix &q, Matrix &p, unsigned long long &computationAmount, unsigned long long &heapOperationAmount1, unsigned long long &heapOperationAmount2) {

    computationAmount = 0;
    heapOperationAmount1 = 0;
    heapOperationAmount2 = 0;

    vector<VectorElement> pNorms(p.rowNum);
    Calculator::calNorms(p, pNorms);
    sort(pNorms.begin(), pNorms.end(), greater<VectorElement>());

    ExtendMatrix<ExtendMatrixRow> preprocessedP;
    preprocessedP.initExtendMatrix(p, pNorms);


    vector<vector<VectorElement> > results(q.rowNum, vector<VectorElement>());
    double qNorm;

    for (int qID = 0; qID < q.rowNum; qID++) {

        const double *qRow = q.getRowPtr(qID);

        vector<VectorElement> &heap = results[qID];
        heap.resize(k);

        for (int pIndex = 0; pIndex < k; pIndex++) {
            const ExtendMatrixRow *pRow = preprocessedP.getRowPtr(pIndex);
            double innerProduct = Calculator::innerProduct(qRow, pRow->rawData, q.colNum);
            heap[pIndex] = VectorElement(pRow->gRowID, innerProduct);
            heapOperationAmount1++;
        }

        make_heap(heap.begin(), heap.end(), greater<VectorElement>());
        double lowerBound = heap.front().data;

        Calculator::calSingleNorm(q.getRowPtr(qID), q.colNum, qNorm);

        for (int pIndex = k; pIndex < preprocessedP.rowNum; pIndex++) {

            const ExtendMatrixRow *pRow = preprocessedP.getRowPtr(pIndex);

            if (pRow->norm * qNorm <= lowerBound) {
                computationAmount += pIndex;
                break;
            }
            else {

                double innerProduct = Calculator::innerProduct(qRow, pRow->rawData, q.colNum);

                if (innerProduct > lowerBound) {
                    pop_heap(heap.begin(), heap.end(), greater<VectorElement>());
                    heap.pop_back();
                    heap.push_back(VectorElement(pRow->gRowID, innerProduct));
                    push_heap(heap.begin(), heap.end(), greater<VectorElement>());
                    lowerBound = heap.front().data;
                    heapOperationAmount2++;
                }
            }
        }
    }

}

void basicNormPruneOnlyScan(const int k, Matrix &q, Matrix &p, unsigned long long &computationAmount, unsigned long long &heapOperationAmount1, unsigned long long &heapOperationAmount2) {

    Monitor t;
    t.start();

    vector<VectorElement> pNorms(p.rowNum);
    Calculator::calNormsForTest(p, pNorms);
    sort(pNorms.begin(), pNorms.end(), greater<VectorElement>());

    ExtendMatrix<ExtendMatrixRow> preprocessedP;
    preprocessedP.initExtendMatrix(p, pNorms);

    vector<vector<VectorElement> > results(q.rowNum, vector<VectorElement>(k));
    volatile double qNorm;
    vector<VectorElement> heap(k);
    volatile double lowerBound = 0;

    for (int i = 0; i < q.rowNum; i++) {
        Calculator::calSingleNormForTest(q.getRowPtr(0), q.colNum, qNorm);
        make_heap(heap.begin(), heap.end(), greater<VectorElement>());
        lowerBound = heap.front().data;
    }

    volatile double innerProduct;

    for (unsigned long long i = 0; i < computationAmount; i++) {
        volatile const double *qRow = q.getRowPtr(0);
        volatile const ExtendMatrixRow *pRow = preprocessedP.getRowPtr(0);
        innerProduct = Calculator::innerProductForTest(qRow, pRow->rawData, q.colNum);
    }

    volatile VectorElement tmp;
    for (unsigned long long i = 0; i < heapOperationAmount1; i++) {
        tmp.id = 0;
        tmp.data = innerProduct;
    }

    for (unsigned long long i = 0; i < heapOperationAmount2; i++) {
        pop_heap(heap.begin(), heap.end(), greater<VectorElement>());
        heap.pop_back();
        heap.push_back(VectorElement(0, innerProduct));
        push_heap(heap.begin(), heap.end(), greater<VectorElement>());
        lowerBound = heap.front().data;
    }

    t.stop();

    cout << "Scan Time: " << t.getElapsedTime() << " secs" << endl;
    cout << lowerBound << endl;
    cout << results[0].size() << endl;
    cout << heap[0].data << endl;

}

void basicNormPruneOnlyCompute(const int qSize, const int dim, int *flag) {

    // Compute how many times we need to compute inner product
    unsigned long long totalIP = 0;
    for (int i = 0; i < qSize; i++) {
//		cout << flag[i] << endl;
        totalIP += flag[i];
    }

    cout << "Computation Amount: " << totalIP << endl;

    double qNorm = 1.0;
    double pNorm = 1.0;
    double threshold = 0.0;
    double *qRowRawData = new double[dim]();
    double *pRowRawData = new double[dim]();

//    srand(time(0));
    for (int i = 0; i < dim; i++) {
//        qRowRawData[i] = rand();
//        pRowRawData[i] = rand();

        qRowRawData[i] = 1.0;
        pRowRawData[i] = 1.0;
    }

    Monitor t;
    t.start();

    // http://stackoverflow.com/questions/15309993/how-do-i-force-the-compiler-not-to-skip-my-function-calls
    // I have test the following for once and the time is 0 sec, while the loop needs more than 0 secs.
    // So the compiler will not ignore the loop and only use the last iteration.
//    volatile double tmp = qNorm * pNorm;
//    if (tmp < threshold) {
//        return;
//    }

//    for (unsigned long long i = 0; i < totalIP; i++) {
//        volatile double tmp = qNorm * pNorm;
//        if (tmp < threshold) {
//            break;
//        }
//    }

    for (unsigned long long i = 0; i < totalIP; i++) {

        volatile double tmp = 1.0 * 0.0;
        if (tmp < 0.0) {
            return;
        }

        volatile double tmp2 = 0;
        for (int j = 0; j < dim; j++) {
//            tmp2 += qRowRawData[j] * pRowRawData[j];
            tmp2 += 1.0 * 0.0;
        }

        if (tmp2 < 0.0) {
            return;
        }

    }

    t.stop();
    cout << "Compute Time: " << t.getElapsedTime() << " secs" << endl;

}

#endif //BASICNORMPRUNE_H