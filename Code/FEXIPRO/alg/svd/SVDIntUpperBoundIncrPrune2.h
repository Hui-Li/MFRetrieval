#ifndef SIRPRUNE_H
#define SIRPRUNE_H

#include "../../util/Base.h"
#include "../../structs/Matrix.h"
#include "../../util/Conf.h"
#include "../../util/SVDUtil.h"
#include "../../util/Monitor.h"
#include "../../util/Calculator.h"
#include "../../structs/SVDIntMatrixRow.h"
#include "../../structs/ExtendMatrix.h"
#include "../../structs/FastHeap.h"
#include "../../structs/SIRMatrixRow.h"

// FEIPR-SIR
class SIRPrune {

private:
    Monitor tt;
    double offlineTime;

    ExtendMatrix<SIRMatrixRow> *preprocessedP;
    Matrix *q;
    Matrix *u;
    int k;
    int checkDim;
    int checkDim2;
    int numOfDim;
    int scalingValue;
    vector<double> addend;
    double minValue;

    // log
#ifdef TIME_IT
    uint64_t counter1 = 0;
    uint64_t counter2 = 0;
    uint64_t counter3 = 0;
    uint64_t counter4 = 0;
    uint64_t counter5 = 0;
#endif

    void transferQ(double &qNorm, double &newSVDQNorm, double &subQNorm, double &subTransformedQNorm,
                   double &sumOfCoordinate, double &leftPartialSumOfCoordinate, const double *qPtr,
                   double *newQ,
                   int *qIntPtr, int &qSumOfCoordinate1, int &qSumOfCoordinate2, double &qRatio1,
                   double &qRatio2);
    void refine(VectorElement *heap, const double qNorm, double *newQ,
                const double subQNorm, const int *newQIntPtr, int qSumOfCoordinate1, int qSumOfCoordinate2,
                double ratio1, double ratio2, const double newSVDQNorm, const double subTransformQNorm,
                const double sumOfQCoordinate, const double leftPartialQSumOfCoordinate);

public:
    SIRPrune(const int k, const int scalingValue, const double SIGMA, Matrix *q, Matrix *p);
    ~SIRPrune();
    void topK();

};

inline void SIRPrune::transferQ(double &qNorm, double &newSVDQNorm, double &subQNorm, double &subTransformedQNorm,
                                double &sumOfCoordinate, double &leftPartialSumOfCoordinate, const double *qPtr,
                                double *newQ,
                                int *qIntPtr, int &qSumOfCoordinate1, int &qSumOfCoordinate2, double &qRatio1,
                                double &qRatio2) {

    double minValue = DBL_MAX;
    double maxValue = -DBL_MAX;
    qNorm = 0;

    for (int colIndex = 0; colIndex < u->colNum; colIndex++) {
        qNorm += qPtr[colIndex] * qPtr[colIndex];
    }
    qNorm = sqrt(qNorm);

    newSVDQNorm = 0;
    for (int rowIndex = u->rowNum - 1; rowIndex >= checkDim; rowIndex--) {

        newQ[rowIndex] = 0;
        const double *uPtr = u->getRowPtr(rowIndex);

        for (int colIndex = 0; colIndex < u->colNum; colIndex++) {
            newQ[rowIndex] += uPtr[colIndex] * qPtr[colIndex];
        }

        if (newQ[rowIndex] < minValue) {
            minValue = newQ[rowIndex];
        }

        if (newQ[rowIndex] > maxValue) {
            maxValue = newQ[rowIndex];
        }

        newSVDQNorm += newQ[rowIndex] * newQ[rowIndex];

    }

    subQNorm = sqrt(newSVDQNorm);

    double absMin = fabs(minValue);
    double absMax = fabs(maxValue); // maxValue can be negative
    double denominator = absMin > absMax ? absMin : absMax;

    if (denominator == 0) {
        qRatio2 = 0;
    } else {
        qRatio2 = scalingValue / denominator;
    }

    qSumOfCoordinate2 = 0; //sumOfCoordinate

    for (int colIndex = q->colNum - 1; colIndex >= checkDim; colIndex--) {
        qIntPtr[colIndex] = floor(newQ[colIndex] * qRatio2);
        qSumOfCoordinate2 += fabs(qIntPtr[colIndex]);
    }

    minValue = DBL_MAX;
    maxValue = -DBL_MAX;

    for (int rowIndex = checkDim - 1; rowIndex >= 0; rowIndex--) {

        newQ[rowIndex] = 0;
        const double *uPtr = u->getRowPtr(rowIndex);

        for (int colIndex = 0; colIndex < u->colNum; colIndex++) {
            newQ[rowIndex] += uPtr[colIndex] * qPtr[colIndex];
        }

        newSVDQNorm += newQ[rowIndex] * newQ[rowIndex];

        if (newQ[rowIndex] < minValue) {
            minValue = newQ[rowIndex];
        }

        if (newQ[rowIndex] > maxValue) {
            maxValue = newQ[rowIndex];
        }

    }

    newSVDQNorm = sqrt(newSVDQNorm);

    subTransformedQNorm = 0;
    sumOfCoordinate = 0;

    double tmp;
    for (int colIndex = numOfDim; colIndex >= checkDim; colIndex--) {
        sumOfCoordinate += addend[colIndex] * newQ[colIndex] / newSVDQNorm;
        tmp = (newQ[colIndex] / newSVDQNorm + addend[colIndex]) * 2;
        subTransformedQNorm += tmp * tmp;
    }
    subTransformedQNorm = sqrt(subTransformedQNorm);
    leftPartialSumOfCoordinate = 0;

    absMin = fabs(minValue);
    absMax = fabs(maxValue); // maxValue can be negative
    denominator = absMin > absMax ? absMin : absMax;

    if (denominator == 0) {
        qRatio1 = 0;
    } else {
        qRatio1 = scalingValue / denominator;
    }

    qSumOfCoordinate1 = 0; //sumOfCoordinate

    for (int colIndex = checkDim - 1; colIndex >= 0; colIndex--) {
        qIntPtr[colIndex] = floor(newQ[colIndex] * qRatio1);
        qSumOfCoordinate1 += fabs(qIntPtr[colIndex]);

        tmp = addend[colIndex] * newQ[colIndex] / newSVDQNorm;
        sumOfCoordinate += tmp;
        leftPartialSumOfCoordinate += tmp;
    }

    qRatio1 = 1 / qRatio1;
    qRatio2 = 1 / qRatio2;

}

inline void SIRPrune::refine(VectorElement *heap, const double qNorm, double *newQ,
                             const double subQNorm, const int *newQIntPtr, int qSumOfCoordinate1, int qSumOfCoordinate2,
                             double ratio1, double ratio2, const double newSVDQNorm, const double subTransformQNorm,
                             const double sumOfQCoordinate, const double leftPartialQSumOfCoordinate) {

    int heapCount = 0;
    for (int rowIndex = 0; rowIndex < k; rowIndex++) {
        const SIRMatrixRow *pRowPtr = preprocessedP->getRowPtr(rowIndex);
        heap_enqueue(Calculator::innerProduct(newQ, pRowPtr->rawData, q->colNum), pRowPtr->gRowID, heap, &heapCount);
    }

    double originalLowerBound = heap[0].data;

    int gRowIDForBoundItem = heap[0].id;
    const SIRMatrixRow *pRowPtr = preprocessedP->getRowPtr(gRowIDForBoundItem);

    double lowerBoundInNewSpace = (originalLowerBound/newSVDQNorm + sumOfQCoordinate + pRowPtr->sumOfCoordinate) * 2 + pRowPtr->partialSumOfCoordinate;

    for (int rowIndex = k; rowIndex < preprocessedP->rowNum; rowIndex++) {

        const SIRMatrixRow *pRowPtr = preprocessedP->getRowPtr(rowIndex);

        if (pRowPtr->norm * qNorm <= originalLowerBound) {
            break;
        }
        else {
#ifdef TIME_IT
            counter1++;
#endif

            int bound = pRowPtr->sumOfCoordinate1 + qSumOfCoordinate1;
            const int* pIntPtr = pRowPtr->iRawData;

            for(int dim = 0; dim < checkDim; dim++) {
                bound += pIntPtr[dim] * newQIntPtr[dim];
            }

            double subNormBound = subQNorm * pRowPtr->subNorm;
            double leftInt = bound * ratio1;
            if(leftInt + subNormBound <= originalLowerBound){
                continue;
            }

#ifdef TIME_IT
            counter2++;
#endif

            bound = pRowPtr->sumOfCoordinate2 + qSumOfCoordinate2;
            for(int dim = checkDim; dim < q->colNum; dim++) {
                bound += pIntPtr[dim] * newQIntPtr[dim];
            }

            if(leftInt + bound * ratio2 <= originalLowerBound) {
                continue;
            }

#ifdef TIME_IT
            counter3++;
#endif

            double innerProduct = 0;
            const double *pPtr = pRowPtr->rawData;
            for(int dim = 0; dim < checkDim; dim++) {
                innerProduct += pPtr[dim] * newQ[dim];
            }

            if(innerProduct + subNormBound <= originalLowerBound){
                continue;
            }

#ifdef TIME_IT
            counter4++;
#endif

            double sublLowerBoundInNewSpace = (innerProduct/newSVDQNorm + leftPartialQSumOfCoordinate + pRowPtr->leftPartialSumOfCoordinate) * 2 + pRowPtr->partialSumOfCoordinate;

            if (sublLowerBoundInNewSpace + pRowPtr->subTransformedSubVNorm * subTransformQNorm <= lowerBoundInNewSpace) {
                continue;
            }

#ifdef TIME_IT
            counter5++;
#endif

            for(int dim = checkDim; dim < q->colNum; dim++) {
                innerProduct += pPtr[dim] * newQ[dim];
            }

            if (innerProduct > originalLowerBound) {

                heap_dequeue(heap, &heapCount);
                heap_enqueue(innerProduct, pRowPtr->gRowID, heap, &heapCount);
                originalLowerBound = heap[0].data;

                int gRowIDForBoundItem = heap[0].id;
                const SIRMatrixRow *pRowPtr = preprocessedP->getRowPtr(gRowIDForBoundItem);
                lowerBoundInNewSpace = (originalLowerBound/newSVDQNorm + sumOfQCoordinate + pRowPtr->sumOfCoordinate) * 2 + pRowPtr->partialSumOfCoordinate;

            }
        }
    }
}

inline SIRPrune::SIRPrune(const int k, const int scalingValue, const double SIGMA, Matrix *q, Matrix *p) {

    mat P_t;
    P_t.load(Conf::pDataPath, csv_ascii);
    P_t = P_t.t();

    tt.start();

    this->q = q;
    this->u = new Matrix();
    this->k = k;
    this->scalingValue = scalingValue;
    this->numOfDim = q->colNum;

    Matrix *v = new Matrix();
//    this->checkDim = SVD(P_t, p->colNum, p->rowNum, *u, *v, SIGMA);

    // add the two extended dimensions
    this->checkDim = SVD(P_t, p->colNum, p->rowNum, *u, *v, addend, SIGMA);
    this->checkDim2 = checkDim + 2;

    preprocessedP = new ExtendMatrix<SIRMatrixRow>();
    preprocessedP->initSIRMatrix(*p, *v, checkDim, checkDim2, scalingValue, addend, minValue);

    delete v;

    tt.stop();
    offlineTime = tt.getElapsedTime();
}

inline SIRPrune::~SIRPrune() {

    if (u) {
        delete u;
    }

    if (preprocessedP) {
        delete preprocessedP;
    }

}

inline void SIRPrune::topK() {

    tt.start();

    VectorElement *results = new VectorElement[q->rowNum * k];

    double *newQ = new double[q->colNum];
    int *newQIntPtr = new int[q->colNum];
    int qSumOfCoordinate1 = 0;
    int qSumOfCoordinate2 = 0;

    double subQNorm = 0;
    double transformSubQNorm = 0;
    double newSVDQNorm = 0;
    double qNorm = 0;
    double sumOfQCoordinate = 0;
    double leftPartialQSumOfCoordinate = 0;

    double qRatio1;
    double qRatio2;
    double ratio1;
    double ratio2;

    const double pRatio1 = preprocessedP->ratio1;
    const double pRatio2 = preprocessedP->ratio2;

    for (int qID = 0; qID < q->rowNum; qID++) {

        VectorElement *heap = &results[qID * k];

        const double *qPtr = q->getRowPtr(qID);

        transferQ(qNorm, newSVDQNorm, subQNorm, transformSubQNorm, sumOfQCoordinate, leftPartialQSumOfCoordinate, qPtr,
                  newQ, newQIntPtr, qSumOfCoordinate1, qSumOfCoordinate2, qRatio1, qRatio2);

        ratio1 = qRatio1 * pRatio1;
        ratio2 = qRatio2 * pRatio2;

        refine(heap, qNorm, newQ, subQNorm, newQIntPtr, qSumOfCoordinate1, qSumOfCoordinate2, ratio1, ratio2,
               newSVDQNorm, transformSubQNorm, sumOfQCoordinate, leftPartialQSumOfCoordinate);

    }

    tt.stop();

    Logger::Log("preprocess time: " + to_string(offlineTime) + " secs");
    Logger::Log("online time: " + to_string(tt.getElapsedTime()) + " secs");

    if (Conf::outputResult) {
        string resultFileName = Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string
                                                                                                              (Conf::k) + "-" + to_string(Conf::scalingValue) + "-" +
                                to_string(Conf::SIGMA) + ".txt";
        FileUtil::outputResult(q->rowNum, k, results, resultFileName);
    }

    delete[] newQ;
    delete[] newQIntPtr;
    delete[] results;
}

#endif //SIRPRUNE_H
