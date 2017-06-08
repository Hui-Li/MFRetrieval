#ifndef SIMDINTUPPERBOUND_H
#define SIMDINTUPPERBOUND_H

#include "../../util/Base.h"
#include "../../structs/Matrix.h"
#include "../../util/Conf.h"
#include "../../util/Monitor.h"
#include "../../util/Calculator.h"
#include "../../util/Logger.h"
#include "../../util/SIMDUtil.h"
#include "../../structs/SIMDIntMatrixRow.h"
#include "../../structs/ExtendMatrix.h"

// SIMD FEIPR-I
// ToDo: It is not working well so far. Using Int8 is not much very fast, and even worse than using int???
class SIMDIntUpperBound {

private:
    Monitor tt;
    double offlineTime;
    Matrix *q;
    ExtendMatrix<SIMDIntMatrixRow> *preprocessedP;
    int k;
    int dimension;
    int scalingValue;
    int lowerDimensionForInt;
    int lowerDimensionForDouble;

#ifdef TIME_IT
    uint64_t counter1;
	uint64_t counter2;
#endif

    double transferQ(const double *qPtr, int16_t *qIntPtr, double &qNorm, int &qSumOfCoordinate);

public:
    SIMDIntUpperBound (const int k, Matrix *q, Matrix *p);
    ~SIMDIntUpperBound();
    void topK();
};


inline double SIMDIntUpperBound::transferQ(const double *qPtr, int16_t *qIntPtr, double &qNorm, int &qSumOfCoordinate){

    double minValue;
    double maxValue;

    SIMDUtil::minMaxValue(dimension, lowerDimensionForDouble, qPtr, maxValue, minValue);

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



inline SIMDIntUpperBound::SIMDIntUpperBound (const int k, Matrix *q, Matrix *p) {

#ifdef __AVX2__
    Logger::Log("Use AVX2: true");
#else
    Logger::Log("Use AVX2: false");
#endif

#ifdef TIME_IT
    counter1 = 0;
    counter2 = 0;
#endif

    tt.start();

    this->q = q;
    this->dimension = p->colNum;
    this->preprocessedP = new ExtendMatrix<SIMDIntMatrixRow>();
    this->k = k;
    this->scalingValue = 127;
    this->lowerDimensionForInt = k / 16 * 16;
    this->lowerDimensionForDouble = k / 4 * 4;

    vector<VectorElement> pNorms(p->rowNum);
    SIMDUtil::doubleNorms(dimension, lowerDimensionForDouble, *p, pNorms);
    sort(pNorms.begin(), pNorms.end(), greater<VectorElement>());

    preprocessedP->initIntExtendMatrixSIMD(*p, pNorms, scalingValue);

    tt.stop();
    offlineTime = tt.getElapsedTime();

}

SIMDIntUpperBound::~SIMDIntUpperBound() {

    if (preprocessedP) {
        delete preprocessedP;
    }

}

inline void SIMDIntUpperBound::topK() {

    tt.start();

    vector<vector<VectorElement> > results(q->rowNum, vector<VectorElement>());

    int16_t *qIntPtr = new int16_t[dimension];
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
            const SIMDIntMatrixRow *pRowPtr = preprocessedP->getRowPtr(pIndex);
            heap[pIndex] = VectorElement(pRowPtr->gRowID, SIMDUtil::doubleDotProduct(dimension, lowerDimensionForDouble, qPtr, pRowPtr->rawData));
        }

        make_heap(heap.begin(), heap.end(), greater<VectorElement>());
        double originalLowerBound = heap.front().data;
        double scaledLowerBound = originalLowerBound * ratio - qSumOfCoordinate;

        for (int pIndex = k; pIndex < preprocessedP->rowNum; pIndex++) {

            const SIMDIntMatrixRow *pRowPtr = preprocessedP->getRowPtr(pIndex);

            if (pRowPtr->norm * qNorm <= originalLowerBound) {
                break;
            } else {

#ifdef TIME_IT
                counter1++;
#endif

                const int16_t *pIntPtr = pRowPtr->iRawData;

                // a*b + a + b + 1
                int bound = SIMDUtil::int32ForIntDotProduct(dimension, lowerDimensionForInt, qIntPtr, pIntPtr) + pRowPtr->sumOfCoordinate;

                if (bound <= scaledLowerBound) {
                    continue;
                }

#ifdef TIME_IT
                counter2++;
#endif

                const double *pPtr = pRowPtr->rawData;
                double innerProduct = SIMDUtil::doubleDotProduct(dimension, lowerDimensionForDouble, qPtr, pRowPtr->rawData);

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
    Logger::Log("Avg Candidate Size after Norm Prune: " + to_string((SIMDIntUpperBound::counter1 + 0.0) / q->rowNum));
    Logger::Log("Avg Candidate Size after Using Int Upper Bound: " + to_string((SIMDIntUpperBound::counter2 + 0.0) / q->rowNum));
#endif

    if (Conf::outputResult) {
        string resultFileName = Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + "-" + to_string(Conf::scalingValue) + ".txt";
        FileUtil::outputResult(k, results, resultFileName);
    }
}

#endif //SIMDINTUPPERBOUND_H
