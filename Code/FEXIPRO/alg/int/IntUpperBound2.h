#ifndef INTUPPERBOUND2_H
#define INTUPPERBOUND2_H

#include "../../util/Base.h"
#include "../../structs/Matrix.h"
#include "../../util/Conf.h"
#include "../../util/Monitor.h"
#include "../../util/Calculator.h"
#include "../../util/Logger.h"
#include "../../structs/IntMatrixRow.h"
#include "../../structs/ExtendMatrix.h"

// FEIPR-I
class IntUpperBound2 {

private:
    Monitor tt;
    double offlineTime;
    Matrix *q;
    Matrix *p;
    ExtendMatrix<IntMatrixRow> *preprocessedP;
    int k;
    int dimension;
    int checkDim;
    int scalingValue;

#ifdef TIME_IT
    uint64_t counter1;
	uint64_t counter2;
    uint64_t counter3;
	uint64_t counter4;
#endif

    void transferQ(const double *qPtr, double *scaledQPtr, int *qIntPtr, double &qNorm,
                                   double &qSubNorm, int &qSumOfCoordinateLeft, int &qSumOfCoordinateRight) ;

public:
    IntUpperBound2(const int k, const int checkDim, const int scalingValue, Matrix *q, Matrix *p);
    ~IntUpperBound2();
    void topK();
};


inline void IntUpperBound2::transferQ(const double *qPtr, double *scaledQPtr, int *qIntPtr, double &qNorm,
                                      double &qSubNorm, int &qSumOfCoordinateLeft, int &qSumOfCoordinateRight) {

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
    qSumOfCoordinateRight = 0; //sumOfCoordinate

    double tmp;

    for (int colIndex = q->colNum-1; colIndex >= checkDim; colIndex--) {
        tmp = qPtr[colIndex] * ratio;
        scaledQPtr[colIndex] = tmp;
        qNorm += scaledQPtr[colIndex] * scaledQPtr[colIndex];
        tmp = floor(tmp);
        qIntPtr[colIndex] = tmp;
        qSumOfCoordinateRight += fabs(tmp);
    }
    qSubNorm = sqrt(qNorm);

    qSumOfCoordinateLeft = 0;
    for (int colIndex = checkDim - 1; colIndex >= 0; colIndex--) {
        tmp = qPtr[colIndex] * ratio;
        scaledQPtr[colIndex] = tmp;
        qNorm += scaledQPtr[colIndex] * scaledQPtr[colIndex];
        tmp = floor(tmp);
        qIntPtr[colIndex] = tmp;
        qSumOfCoordinateLeft += fabs(tmp);
    }
    qNorm = sqrt(qNorm);

}



inline IntUpperBound2::IntUpperBound2(const int k, const int checkDim, const int scalingValue, Matrix *q, Matrix *p) {

#ifdef TIME_IT
    counter1 = 0;
    counter2 = 0;
    counter3 = 0;
    counter4 = 0;
#endif

    tt.start();

    this->q = q;
    this->p = p;
    this->dimension = p->colNum;
    this->checkDim = checkDim;
    this->preprocessedP = new ExtendMatrix<IntMatrixRow>();
    this->k = k;
    this->scalingValue = scalingValue;

    vector<VectorElement> pNorms(p->rowNum);
    Calculator::calNorms(*p, pNorms);
    sort(pNorms.begin(), pNorms.end(), greater<VectorElement>());
    preprocessedP->initIntExtendMatrix2(*p, pNorms, scalingValue, checkDim);

    tt.stop();
    offlineTime = tt.getElapsedTime();

}

IntUpperBound2::~IntUpperBound2() {

    if (preprocessedP) {
        delete preprocessedP;
    }

}

inline void IntUpperBound2::topK() {

    tt.start();

    vector<vector<VectorElement> > results(q->rowNum, vector<VectorElement>());

    int *qIntPtr = new int[dimension];

    double qNorm;
    double qSubNorm;
    int qSumOfCoordinateLeft;
    int qSumOfCoordinateRight;
    double *scaledQPtr = new double[this->dimension];

    for (int qID = 0; qID < q->rowNum; qID++) {

        // Transfer q to int vector
        const double *qPtr = q->getRowPtr(qID);

        transferQ(qPtr, scaledQPtr, qIntPtr, qNorm, qSubNorm, qSumOfCoordinateLeft, qSumOfCoordinateRight);

        vector<VectorElement> &heap = results[qID];
        heap.resize(k);

        for (int pIndex = 0; pIndex < k; pIndex++) {
            const IntMatrixRow *pRowPtr = preprocessedP->getRowPtr(pIndex);
            heap[pIndex] = VectorElement(pRowPtr->gRowID, Calculator::innerProduct(scaledQPtr, pRowPtr->rawData, q->colNum));
        }

        make_heap(heap.begin(), heap.end(), greater<VectorElement>());
        double lowerBound = heap.front().data;

        for (int pIndex = k; pIndex < preprocessedP->rowNum; pIndex++) {

            const IntMatrixRow *pRowPtr = preprocessedP->getRowPtr(pIndex);

            if (pRowPtr->norm * qNorm <= lowerBound) {
                break;
            } else {

#ifdef TIME_IT
                counter1++;
#endif

                const int *pIntPtr = pRowPtr->iRawData;
                int intBound = pRowPtr->sumOfCoordinateLeft + qSumOfCoordinateLeft;
                // a*b + a + b + 1

                for (int colIndex = 0; colIndex < checkDim; colIndex++) {
                    intBound += qIntPtr[colIndex] * pIntPtr[colIndex];
                }

                double normBound = pRowPtr->subNorm * qSubNorm;
                if (intBound + normBound <= lowerBound) {
                    continue;
                }

#ifdef TIME_IT
                counter2++;
#endif

                for (int colIndex = checkDim; colIndex < dimension; colIndex++) {
                    intBound += qIntPtr[colIndex] * pIntPtr[colIndex];
                }

                intBound += pRowPtr->sumOfCoordinateRight + qSumOfCoordinateRight;
                if (intBound <= lowerBound) {
                    continue;
                }

#ifdef TIME_IT
                counter3++;
#endif

                const double *pScaledPtr = pRowPtr->rawData;
                double innerProduct = 0;
                for (int colIndex = 0; colIndex < checkDim; colIndex++) {
                    innerProduct += scaledQPtr[colIndex] * pScaledPtr[colIndex];
                }

                if (innerProduct + normBound <= lowerBound){
                    continue;
                }

#ifdef TIME_IT
                counter4++;
#endif

                for (int colIndex = checkDim; colIndex < dimension; colIndex++) {
                    innerProduct += scaledQPtr[colIndex] * pScaledPtr[colIndex];
                }

                if (innerProduct > lowerBound) {
                    pop_heap(heap.begin(), heap.end(), greater<VectorElement>());
                    heap.pop_back();
                    heap.push_back(VectorElement(pRowPtr->gRowID, innerProduct));
                    push_heap(heap.begin(), heap.end(), greater<VectorElement>());
                    lowerBound = heap.front().data;
                }

            }

        }
    }

    tt.stop();

    delete[] scaledQPtr;
    delete[] qIntPtr;

    Logger::Log("preprocess time: " + to_string(offlineTime) + " secs");
    Logger::Log("online time: " + to_string(tt.getElapsedTime()) + " secs");

#ifdef TIME_IT
    Logger::Log("counter1: " + to_string((IntUpperBound2::counter1 + 0.0) / q->rowNum));
    Logger::Log("counter2: " + to_string((IntUpperBound2::counter2 + 0.0) / q->rowNum));
    Logger::Log("counter3: " + to_string((IntUpperBound2::counter3 + 0.0) / q->rowNum));
    Logger::Log("counter4: " + to_string((IntUpperBound2::counter4 + 0.0) / q->rowNum));
#endif

    if (Conf::outputResult) {
        string resultFileName = Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + "-" + to_string(Conf::scalingValue) + ".txt";
        vector<vector<VectorElement> > originalResults(q->rowNum, vector<VectorElement>(k));

        for (int qIndex = 0; qIndex < results.size(); qIndex++) {
            for (int i = 0; i < k; i++) {
                originalResults[qIndex][i] = VectorElement(results[qIndex][i].id, Calculator::innerProduct(q->getRowPtr(qIndex), p->getRowPtr(results[qIndex][i].id), dimension));
            }
        }
        FileUtil::outputResult(k, originalResults, resultFileName);
    }

}

#endif //INTUPPERBOUND2_H
