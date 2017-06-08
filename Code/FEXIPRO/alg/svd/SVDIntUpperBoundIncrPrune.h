#ifndef SVDINTUPPERBOUNDINCRPRUNE_H
#define SVDINTUPPERBOUNDINCRPRUNE_H

#include "../../util/Base.h"
#include "../../structs/Matrix.h"
#include "../../util/Conf.h"
#include "../../util/SVDUtil.h"
#include "../../util/Monitor.h"
#include "../../util/Calculator.h"
#include "../../structs/SVDIntMatrixRow.h"
#include "../../structs/ExtendMatrix.h"
#include "../../structs/FastHeap.h"

// FEIPR-SI
class SVDIntUpperBoundIncrPrune {

private:
    Monitor tt;
    double offlineTime;

    ExtendMatrix<SVDIntMatrixRow> *preprocessedP;
    Matrix *q;
    Matrix *u;
    int k;
    int checkDim;
    int scalingValue;

    // log
#ifdef TIME_IT
    uint64_t counter1 = 0;
    uint64_t counter2 = 0;
    uint64_t counter3 = 0;
    uint64_t counter4 = 0;
#endif

    void transferQ(double &qNorm, double &subQNorm, const double *qPtr, double *newQ,
                     int *qIntPtr, int &qSumOfCoordinate1, int &qSumOfCoordinate2, double &qRatio1, double &qRatio2);
    void refine(VectorElement *heap, const double qNorm, double *newQ,
                const double subQNorm, const int *newQIntPtr, int qSumOfCoordinate1, int qSumOfCoordinate2, double ratio1, double ratio2);

public:
    SVDIntUpperBoundIncrPrune(const int k, const int scalingValue, const double SIGMA, Matrix *q, Matrix *p);
    ~SVDIntUpperBoundIncrPrune();
    void topK();

};

inline void SVDIntUpperBoundIncrPrune::transferQ(double &qNorm, double &subQNorm, const double *qPtr, double *newQ,
                        int *qIntPtr, int &qSumOfCoordinate1, int &qSumOfCoordinate2, double &qRatio1, double &qRatio2) {

    double minValue = DBL_MAX;
    double maxValue = -DBL_MAX;
    qNorm = 0;

    for (int colIndex = 0; colIndex < u->colNum; colIndex++) {
        qNorm += qPtr[colIndex] * qPtr[colIndex];
    }
    qNorm = sqrt(qNorm);

    subQNorm = 0;
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

        subQNorm += newQ[rowIndex] * newQ[rowIndex];

    }

    subQNorm = sqrt(subQNorm);

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

        if (newQ[rowIndex] < minValue) {
            minValue = newQ[rowIndex];
        }

        if (newQ[rowIndex] > maxValue) {
            maxValue = newQ[rowIndex];
        }

    }

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
    }

    qRatio1 = 1 / qRatio1;
    qRatio2 = 1 / qRatio2;

}

inline void SVDIntUpperBoundIncrPrune::refine(VectorElement *heap, const double qNorm, double *newQ,
                   const double subQNorm, const int *newQIntPtr, int qSumOfCoordinate1, int qSumOfCoordinate2, double ratio1, double ratio2) {

    int heapCount = 0;
    for (int rowIndex = 0; rowIndex < k; rowIndex++) {
        const SVDIntMatrixRow *pRowPtr = preprocessedP->getRowPtr(rowIndex);
        heap_enqueue(Calculator::innerProduct(newQ, pRowPtr->rawData, q->colNum), pRowPtr->gRowID, heap, &heapCount);
    }

    double originalLowerBound = heap[0].data;

    for (int rowIndex = k; rowIndex < preprocessedP->rowNum; rowIndex++) {

        const SVDIntMatrixRow *pRowPtr = preprocessedP->getRowPtr(rowIndex);

        if (pRowPtr->norm * qNorm <= originalLowerBound) {
            break;
        }
        else {
#ifdef TIME_IT
            counter1++;
#endif

            // svd incr

            int bound = pRowPtr->sumOfCoordinate1 + qSumOfCoordinate1;
            const int* pIntPtr = pRowPtr->iRawData;

            for(int dim = 0; dim < checkDim; dim++) {
                bound += pIntPtr[dim] * newQIntPtr[dim];
            }

            double subNormBound = subQNorm * pRowPtr->subNorm;
            if(bound * ratio1 + subNormBound <= originalLowerBound){
                continue;
            }

#ifdef TIME_IT
            counter2++;
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
            counter3++;
#endif

            bound = pRowPtr->sumOfCoordinate2 + qSumOfCoordinate2;
            for(int dim = checkDim; dim < q->colNum; dim++) {
                bound += pIntPtr[dim] * newQIntPtr[dim];
            }

            if(innerProduct + bound * ratio2 <= originalLowerBound){
                continue;
            }

#ifdef TIME_IT
            counter4++;
#endif

            // all bounds fail, calculate the exact value.
            for(int dim = checkDim; dim < q->colNum; dim++) {
                innerProduct += pPtr[dim] * newQ[dim];
            }

            if (innerProduct > originalLowerBound) {

                heap_dequeue(heap, &heapCount);
                heap_enqueue(innerProduct, pRowPtr->gRowID, heap, &heapCount);
                originalLowerBound = heap[0].data;
            }
        }
    }
}

inline SVDIntUpperBoundIncrPrune::SVDIntUpperBoundIncrPrune(const int k, const int scalingValue, const double SIGMA, Matrix *q, Matrix *p) {

    mat P_t;
    P_t.load(Conf::pDataPath, csv_ascii);
    P_t = P_t.t();

    tt.start();

    this->q = q;
    this->u = new Matrix();
    this->k = k;
    this->scalingValue = scalingValue;

    Matrix *v = new Matrix();
    this->checkDim = SVD(P_t, p->colNum, p->rowNum, *u, *v, SIGMA);

    vector<VectorElement> pNorms(p->rowNum);
    Calculator::calNorms(*p, pNorms);
    sort(pNorms.begin(), pNorms.end(), greater<VectorElement>());

    preprocessedP = new ExtendMatrix<SVDIntMatrixRow>();
    preprocessedP->initSVDIntExtendMatrix(*v, pNorms, checkDim, scalingValue);

    delete v;

    tt.stop();
    offlineTime = tt.getElapsedTime();
}

inline SVDIntUpperBoundIncrPrune::~SVDIntUpperBoundIncrPrune() {

    if (u) {
        delete u;
    }

    if (preprocessedP) {
        delete preprocessedP;
    }

}

inline void SVDIntUpperBoundIncrPrune::topK() {

    tt.start();

    VectorElement *results = new VectorElement[q->rowNum * k];

    double *newQ = new double[q->colNum];
    int *newQIntPtr = new int[q->colNum];
    int qSumOfCoordinate1 = 0;
    int qSumOfCoordinate2 = 0;

    double subQNorm = 0;
    double qRatio1;
    double qRatio2;
    double ratio1;
    double ratio2;
    double qNorm;
    const double pRatio1 = preprocessedP->ratio1;
    const double pRatio2 = preprocessedP->ratio2;

    for (int qID = 0; qID < q->rowNum; qID++) {

        // step 1: transfer q
        VectorElement *heap = &results[qID * k];

        const double *qPtr = q->getRowPtr(qID);

        transferQ(qNorm, subQNorm, qPtr, newQ, newQIntPtr, qSumOfCoordinate1, qSumOfCoordinate2, qRatio1, qRatio2);
        ratio1 = qRatio1 * pRatio1;
        ratio2 = qRatio2 * pRatio2;

        // step 3: refine
        refine(heap, qNorm, newQ, subQNorm, newQIntPtr, qSumOfCoordinate1, qSumOfCoordinate2, ratio1, ratio2);
//            cout << heap[0].id << "," << heap[0].data << endl;
//            exit(0);
    }

    tt.stop();

#ifdef TIME_IT
    Logger::Log("counter1: " + to_string((counter1 + 0.0) / q->rowNum));
        Logger::Log("counter2: " + to_string((counter2 + 0.0) / q->rowNum));
        Logger::Log("counter3: " + to_string((counter3 + 0.0) / q->rowNum));
        Logger::Log("counter4: " + to_string((counter4 + 0.0) / q->rowNum));
#endif

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

#endif //SVDINTUPPERBOUNDINCRPRUNE_H
