#ifndef MAPINTBASICPRUNE_H
#define MAPINTBASICPRUNE_H

#include "../../util/Base.h"
#include "../../structs/Matrix.h"
#include "../../util/Calculator.h"
#include "../../structs/IntMatrixRow.h"
#include "../../structs/ExtendMatrix.h"

// map to int
class MapIntNormPune {

private:

    Matrix *q;
    ExtendMatrix<IntMatrixRow> *preprocessedP;
    int k;
    int dim;
    int mapMaxValue;

    inline double transferQ(const double *qPtr, int *qIntPtr, double &qNorm, int &qSumOfCoordinate);

public:
    inline MapIntNormPune (const int k, const int mapMaxValue, Matrix *q, Matrix *p);
    ~MapIntNormPune();
    inline void topK();
};


inline double MapIntNormPune::transferQ(const double *qPtr, int *qIntPtr, double &qNorm, int &qSumOfCoordinate){

    double minValue = DBL_MAX;
    double maxValue = -DBL_MAX;

    for (int colIndex = 0; colIndex < dim; colIndex++) {

        if(qPtr[colIndex] < minValue) {
            minValue = qPtr[colIndex];
        }

        if(qPtr[colIndex] > maxValue) {
            maxValue = qPtr[colIndex];
        }
    }

    double absMin = fabs(minValue);
    double absMax = fabs(maxValue); // maxValue can be negative
    double denominator = absMin > absMax ? absMin : absMax;
    double ratio;
    if(denominator==0){
        ratio = 0;
    } else {
        ratio = this->mapMaxValue / denominator;
    }

    qNorm = 0;
    qSumOfCoordinate = 0; //sumOfCoordinate

    for (int colIndex = 0; colIndex < q->colNum; colIndex++) {
        qNorm += qPtr[colIndex] * qPtr[colIndex];
        qIntPtr[colIndex] = floor(qPtr[colIndex] * ratio);
        qSumOfCoordinate += fabs(qIntPtr[colIndex]);
    }

    qNorm = sqrt(qNorm);

    return ratio;
}




inline MapIntNormPune::MapIntNormPune (const int k, const int mapMaxValue, Matrix *q, Matrix *p) {

    this->q = q;
    this->dim = p->colNum;
    this->preprocessedP = new ExtendMatrix<IntMatrixRow>();
    this->k = k;
    this->mapMaxValue = mapMaxValue;

    vector<VectorElement> pNorms(p->rowNum);
    Calculator::calNorms(*p, pNorms);
    sort(pNorms.begin(), pNorms.end(), greater<VectorElement>());
    preprocessedP->initIntExtendMatrix(*p, pNorms, mapMaxValue);

}

MapIntNormPune::~MapIntNormPune() {

    if (preprocessedP) {
        delete preprocessedP;
    }

}


// bound for those cannot be pruned
inline void MapIntNormPune::topK() {

    vector<vector<VectorElement> > results(q->rowNum, vector<VectorElement>());

    int *qIntPtr = new int[dim];
    double qRatio;
    double qNorm;
    int qSumOfCoordinate;
    double pRatio = preprocessedP->ratio;
    double ratio;

//    double sum1 = 0;
    double sum2 = 0;
    uint64_t count = 0;

    for (int qID = 0; qID < q->rowNum; qID++) {

        // Transfer q to int vector
        const double *qPtr = q->getRowPtr(qID);

        qRatio = transferQ(qPtr, qIntPtr, qNorm, qSumOfCoordinate);;
        ratio = qRatio * pRatio;
        vector<VectorElement> &heap = results[qID];
        heap.resize(k);

        for (int pIndex = 0; pIndex < k; pIndex++) {

            const IntMatrixRow *pRowPtr = preprocessedP->getRowPtr(pIndex);

            heap[pIndex] = VectorElement(pRowPtr->gRowID, Calculator::innerProduct(qPtr, pRowPtr->rawData, q->colNum));

        }

        make_heap(heap.begin(), heap.end(), greater<VectorElement>());
        double originalLowerBound = heap.front().data;
        double scaledLowerBound = originalLowerBound * ratio;

        for (int pIndex = k; pIndex < preprocessedP->rowNum; pIndex++) {

            const IntMatrixRow *pRowPtr = preprocessedP->getRowPtr(pIndex);

            if (pRowPtr->norm * qNorm <= originalLowerBound) {
                break;
            } else {

                const int *pIntPtr = pRowPtr->iRawData;
                int bound = qSumOfCoordinate + pRowPtr->sumOfCoordinate;
                // a*b + a + b + 1

                for (int colIndex = 0; colIndex < dim; colIndex++) {
                    bound += qIntPtr[colIndex] * pIntPtr[colIndex];
                }

                const double *pPtr = pRowPtr->rawData;
                double innerProduct = Calculator::innerProduct(qPtr, pPtr, dim);

                double tmp = bound / (innerProduct * ratio);

//                // overflow
//                if (!(isinf(tmp) || tmp < 0)) {
//                    sum1 += tmp;
//                }
//
                if (bound <= scaledLowerBound) {
                    continue;
                }

                count++;

                // overflow
                if (!(std::isinf(tmp) || tmp < 0)) {
                    sum2 += tmp;
                }

                if (innerProduct > originalLowerBound) {
                    pop_heap(heap.begin(), heap.end(), greater<VectorElement>());
                    heap.pop_back();
                    heap.push_back(VectorElement(pRowPtr->gRowID, innerProduct));
                    push_heap(heap.begin(), heap.end(), greater<VectorElement>());
                    originalLowerBound = heap.front().data;
                    scaledLowerBound = originalLowerBound * ratio;
                }

            }

        }
    }

    cout << "avg candidate size: " << (count+0.1) / q->rowNum << endl;
//    cout << "Bound Ratio for all: " << sum1 / q->rowNum << endl;
    cout << "Bound Ratio for those cannot be pruned: " << sum2 / count << endl;

}

#endif //MAPINTBASICPRUNE_H
