#ifndef INTUPPERBOUND_H
#define INTUPPERBOUND_H

#include "../../util/Base.h"
#include "../../structs/Matrix.h"
#include "../../util/Conf.h"
#include "../../util/Monitor.h"
#include "../../util/Calculator.h"
#include "../../util/Logger.h"
#include "../../structs/IntMatrixRow.h"
#include "../../structs/ExtendMatrix.h"

// FEIPR-I
class IntUpperBound {

private:
    Monitor tt;
    double offlineTime;
    Matrix *q;
    ExtendMatrix<IntMatrixRow> *preprocessedP;
    int k;
    int dimension;
    int scalingValue;

#ifdef TIME_IT
    uint64_t counter1;
	uint64_t counter2;
#endif

    double transferQ(const double *qPtr, int *qIntPtr, double &qNorm, int &qSumOfCoordinate);

public:
    IntUpperBound (const int k, const int scalingValue, Matrix *q, Matrix *p);
    ~IntUpperBound();
    void topK();
};


inline double IntUpperBound::transferQ(const double *qPtr, int *qIntPtr, double &qNorm, int &qSumOfCoordinate){

    double minValue = DBL_MAX;
    double maxValue = -DBL_MAX;

    for (int colIndex = 0; colIndex < dimension; colIndex++) {

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
        ratio = this->scalingValue / denominator;
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



inline IntUpperBound::IntUpperBound (const int k, const int scalingValue, Matrix *q, Matrix *p) {

#ifdef TIME_IT
    counter1 = 0;
    counter2 = 0;
#endif

    tt.start();

    this->q = q;
    this->dimension = p->colNum;
    this->preprocessedP = new ExtendMatrix<IntMatrixRow>();
    this->k = k;
    this->scalingValue = scalingValue;

    vector<VectorElement> pNorms(p->rowNum);
    Calculator::calNorms(*p, pNorms);
    sort(pNorms.begin(), pNorms.end(), greater<VectorElement>());
    preprocessedP->initIntExtendMatrix(*p, pNorms, scalingValue);

    tt.stop();
    offlineTime = tt.getElapsedTime();

}

IntUpperBound::~IntUpperBound() {

    if (preprocessedP) {
        delete preprocessedP;
    }

}

inline void IntUpperBound::topK() {

    tt.start();

    vector<vector<VectorElement> > results(q->rowNum, vector<VectorElement>());

    int *qIntPtr = new int[dimension];
    double qRatio;
    double qNorm;
    int qSumOfCoordinate;
    double pRatio = preprocessedP->ratio;
    double ratio;

    for (int qID = 0; qID < q->rowNum; qID++) {

        // Transfer q to int vector
        const double *qPtr = q->getRowPtr(qID);

        qRatio = transferQ(qPtr, qIntPtr, qNorm, qSumOfCoordinate);
        ratio = qRatio * pRatio;
        vector<VectorElement> &heap = results[qID];
        heap.resize(k);

        for (int pIndex = 0; pIndex < k; pIndex++) {
            const IntMatrixRow *pRowPtr = preprocessedP->getRowPtr(pIndex);
            heap[pIndex] = VectorElement(pRowPtr->gRowID, Calculator::innerProduct(qPtr, pRowPtr->rawData, q->colNum));
        }

        make_heap(heap.begin(), heap.end(), greater<VectorElement>());
        double originalLowerBound = heap.front().data;
        double scaledLowerBound = originalLowerBound * ratio - qSumOfCoordinate;

        for (int pIndex = k; pIndex < preprocessedP->rowNum; pIndex++) {

            const IntMatrixRow *pRowPtr = preprocessedP->getRowPtr(pIndex);

            if (pRowPtr->norm * qNorm <= originalLowerBound) {
                break;
            } else {

#ifdef TIME_IT
                counter1++;
#endif

                const int *pIntPtr = pRowPtr->iRawData;
                int bound = 0;
                // a*b + a + b + 1

                for (int colIndex = 0; colIndex < dimension; colIndex++) {
                    bound += qIntPtr[colIndex] * pIntPtr[colIndex];
                }

                bound += pRowPtr->sumOfCoordinate;
                if (bound <= scaledLowerBound) {
                    continue;
                }

#ifdef TIME_IT
                counter2++;
#endif

                const double *pPtr = pRowPtr->rawData;
                double innerProduct = Calculator::innerProduct(qPtr, pPtr, dimension);

                if (innerProduct > originalLowerBound) {
                    pop_heap(heap.begin(), heap.end(), greater<VectorElement>());
                    heap.pop_back();
                    heap.push_back(VectorElement(pRowPtr->gRowID, innerProduct));
                    push_heap(heap.begin(), heap.end(), greater<VectorElement>());
                    originalLowerBound = heap.front().data;
                    scaledLowerBound = originalLowerBound * ratio - qSumOfCoordinate;
                }

            }

        }
    }

    tt.stop();

    delete[] qIntPtr;

    Logger::Log("preprocess time: " + to_string(offlineTime) + " secs");
    Logger::Log("online time: " + to_string(tt.getElapsedTime()) + " secs");

#ifdef TIME_IT
    Logger::Log("Avg Candidate Size after Norm Prune: " + to_string((IntUpperBound::counter1 + 0.0) / q->rowNum));
    Logger::Log("Avg Candidate Size after Using Int Upper Bound: " + to_string((IntUpperBound::counter2 + 0.0) / q->rowNum));
#endif

    if (Conf::outputResult) {
        string resultFileName = Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + "-" + to_string(Conf::scalingValue) + ".txt";
        FileUtil::outputResult(k, results, resultFileName);
    }
}

#endif //INTUPPERBOUND_H
