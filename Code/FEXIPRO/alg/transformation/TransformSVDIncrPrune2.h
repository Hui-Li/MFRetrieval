#ifndef TRANSFORMSVDINCRPRUNE2_H
#define TRANSFORMSVDINCRPRUNE2_H

#include "../../util/Base.h"
#include "../../structs/Matrix.h"
#include "../../util/Conf.h"
#include "../../util/SVDUtil.h"
#include "../../util/Monitor.h"
#include "../../util/Calculator.h"
#include "../../structs/RedSVDMatrixRow2.h"

// SR
// svd incr first done over original svd vector, then transformed vector
class TransformSVDIncrPrune2 {
private:
    Monitor tt;
    double offlineTime;

    ExtendMatrix<RedSVDMatrixRow2> *preprocessedP;
    Matrix *p;
    Matrix *q;
    Matrix *u;
    int k;
    int checkDim;
    int checkDim2;
    int numOfNewDim;
    int numOfDim;
    vector<double> addend;
    double minValue;

    // log
#ifdef TIME_IT
    uint64_t counter1 = 0;
    uint64_t counter2 = 0;
    uint64_t counter3 = 0;
#endif

    void transferQ(double &qNorm, double &newSVDQNorm, double &subQNorm, double &subTransformedQNorm,
                   double &sumOfCoordinate, double &leftPartialSumOfCoordinate, const double *qPtr, double *newQ);
    void refine(vector<VectorElement> &heap,
                double *newQ, const double qNorm,
                const double newSVDQNorm, const double &subQNorm, const double &subTransformQNorm,
                const double sumOfQCoordinate, const double leftPartialQSumOfCoordinate);
public:
    TransformSVDIncrPrune2(const int k, const double SIGMA, Matrix *q, Matrix *p);
    ~TransformSVDIncrPrune2();
    void topK();

};

inline void TransformSVDIncrPrune2::transferQ(double &qNorm, double &newSVDQNorm, double &subQNorm, double &subTransformedQNorm,
                                              double &sumOfCoordinate, double &leftPartialSumOfCoordinate, const double *qPtr, double *newQ) {

    Calculator::calSingleNorm(qPtr, numOfDim, qNorm);

    newSVDQNorm = 0;
    for (int rowIndex = u->rowNum - 1; rowIndex >= checkDim; rowIndex--) {
        newQ[rowIndex] = 0;
        const double *uPtr = u->getRowPtr(rowIndex);

        for (int colIndex = 0; colIndex < u->colNum; colIndex++) {
            newQ[rowIndex] += uPtr[colIndex] * qPtr[colIndex];
        }
        newSVDQNorm += newQ[rowIndex] * newQ[rowIndex];
    }

    subQNorm = sqrt(newSVDQNorm);

    for (int rowIndex = checkDim - 1; rowIndex >= 0; rowIndex--) {
        newQ[rowIndex] = 0;
        const double *uPtr = u->getRowPtr(rowIndex);

        for (int colIndex = 0; colIndex < u->colNum; colIndex++) {
            newQ[rowIndex] += uPtr[colIndex] * qPtr[colIndex];
        }
        newSVDQNorm += newQ[rowIndex] * newQ[rowIndex];
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

    for (int colIndex = checkDim-1; colIndex >= 0; colIndex--) {
        tmp = addend[colIndex] * newQ[colIndex] / newSVDQNorm;
        sumOfCoordinate += tmp;
        leftPartialSumOfCoordinate += tmp;
    }
}

inline void TransformSVDIncrPrune2::refine(vector<VectorElement> &heap,
                                           double *newQ, const double qNorm,
                                           const double newSVDQNorm, const double &subQNorm,
                                           const double &subTransformQNorm,
                                           const double sumOfQCoordinate, const double leftPartialQSumOfCoordinate) {
    
    for (int rowIndex = 0; rowIndex < k; rowIndex++) {
        const RedSVDMatrixRow2 *pRowPtr = preprocessedP->getRowPtr(rowIndex);
        double innerProduct = Calculator::innerProduct(newQ, pRowPtr->rawData, numOfDim);
        heap[rowIndex] = VectorElement(pRowPtr->gRowID, innerProduct);
    }

    make_heap(heap.begin(), heap.end(), greater<VectorElement>());
    double originalLowerBound = heap.front().data;

    int gRowIDForBoundItem = heap.front().id;
    const RedSVDMatrixRow2 *pRowPtr = preprocessedP->getRowPtr(gRowIDForBoundItem);

    double lowerBoundInNewSpace = (originalLowerBound/newSVDQNorm + sumOfQCoordinate + pRowPtr->sumOfCoordinate) * 2 + pRowPtr->partialSumOfCoordinate;
//    lowerBoundInNewSpace = isnan(lowerBoundInNewSpace)?0:lowerBoundInNewSpace;

    for (int rowIndex = k; rowIndex < preprocessedP->rowNum; rowIndex++) {

        const RedSVDMatrixRow2 *pRowPtr = preprocessedP->getRowPtr(rowIndex);

        if (pRowPtr->norm * qNorm <= originalLowerBound) {
            break;
        }
        else {
#ifdef TIME_IT
            counter1++;
#endif

            const double *pPtr = pRowPtr->rawData;

            double innerProduct = 0;
            for (int colIndex = 0; colIndex < checkDim; colIndex++) {
                innerProduct += newQ[colIndex] * pPtr[colIndex];
            }

            if (innerProduct + pRowPtr->subVNorm * subQNorm <= originalLowerBound) {
                continue;
            }


#ifdef TIME_IT
            counter2++;
#endif

            double sublLowerBoundInNewSpace = (innerProduct/newSVDQNorm + leftPartialQSumOfCoordinate + pRowPtr->leftPartialSumOfCoordinate) * 2 + pRowPtr->partialSumOfCoordinate;

            if (sublLowerBoundInNewSpace + pRowPtr->subTransformedSubVNorm * subTransformQNorm <= lowerBoundInNewSpace) {
                continue;
            }

#ifdef TIME_IT
            counter3++;
#endif

            for (int colIndex = checkDim; colIndex < numOfDim; colIndex++) {
                innerProduct += newQ[colIndex] * pPtr[colIndex];
            }

            if (innerProduct > originalLowerBound) {
                pop_heap(heap.begin(), heap.end(), greater<VectorElement>());
                heap.pop_back();
                heap.push_back(VectorElement(pRowPtr->gRowID, innerProduct));
                push_heap(heap.begin(), heap.end(), greater<VectorElement>());

                originalLowerBound = heap.front().data;

                int gRowIDForBoundItem = heap.front().id;
                const RedSVDMatrixRow2 *pRowPtr = preprocessedP->getRowPtr(gRowIDForBoundItem);
                lowerBoundInNewSpace = (originalLowerBound/newSVDQNorm + sumOfQCoordinate + pRowPtr->sumOfCoordinate) * 2 + pRowPtr->partialSumOfCoordinate;
//                lowerBoundInNewSpace = isnan(lowerBoundInNewSpace)?0:lowerBoundInNewSpace;
            }

        }
    }

}


inline TransformSVDIncrPrune2::TransformSVDIncrPrune2(const int k, const double SIGMA, Matrix *q, Matrix *p) {

    tt.start();

    this->q = q;
    this->p = p;
    this->u = new Matrix();
    this->k = k;
    this->numOfDim = q->colNum;
    this->numOfNewDim = q->colNum + 2;

    mat P_t(p->rawData, p->colNum, p->rowNum, false, true);

    Matrix *v = new Matrix();

    // add the two extended dimensions
    this->checkDim = SVD(P_t, p->colNum, p->rowNum, *u, *v, addend, SIGMA);
    this->checkDim2 = checkDim + 2;

    preprocessedP = new ExtendMatrix<RedSVDMatrixRow2>();
    preprocessedP->initSVDTransformExtendMatrix2(*p, *v, checkDim2, addend, minValue);

    delete v;

    tt.stop();
    offlineTime = tt.getElapsedTime();
}

inline TransformSVDIncrPrune2::~TransformSVDIncrPrune2() {

    if (u) {
        delete u;
    }

    if (preprocessedP) {
        delete preprocessedP;
    }

}

inline void TransformSVDIncrPrune2::topK() {

    tt.start();

    vector<vector<VectorElement> > results(q->rowNum, vector<VectorElement>());

    double *newQ = new double[numOfDim];
    double subQNorm = 0;
    double transformSubQNorm = 0;
    double newSVDQNorm = 0;
    double qNorm = 0;
    double sumOfQCoordinate = 0;
    double leftPartialQSumOfCoordinate = 0;

    for (int qID = 0; qID < q->rowNum; qID++) {
        vector<VectorElement> &heap = results[qID];
        heap.resize(k);

        const double *qPtr = q->getRowPtr(qID);

        transferQ(qNorm, newSVDQNorm, subQNorm, transformSubQNorm, sumOfQCoordinate, leftPartialQSumOfCoordinate, qPtr, newQ);

        refine(heap, newQ, qNorm, newSVDQNorm, subQNorm, transformSubQNorm, sumOfQCoordinate, leftPartialQSumOfCoordinate);
    }


    delete[] newQ;
    tt.stop();

#ifdef TIME_IT
    Logger::Log("Avg Num of p which can pass first Cauchy Schwarz inequality check: " + to_string((counter1 + 0.0)/ q
    ->rowNum));
    Logger::Log("Avg Num of p which can pass increamental prune on svd vectors: " + to_string((counter2 + 0.0)/ q
    ->rowNum));
    Logger::Log("Avg Num of p need to be calculated exactly: " + to_string((counter3 + 0.0) / q->rowNum));
#endif

    Logger::Log("preprocess time: " + to_string(offlineTime) + " secs");
    Logger::Log("online time: " + to_string(tt.getElapsedTime()) + " secs");

    if (Conf::outputResult) {
        string resultFileName = Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(k) + ".txt";
        FileUtil::outputResult(k, results, resultFileName);
    }
}

#endif //TRANSFORMSVDINCRPRUNE2_H
