#ifndef TRANSFORMSVDINCRPRUNE_H
#define TRANSFORMSVDINCRPRUNE_H

#include "../../util/Base.h"
#include "../../structs/Matrix.h"
#include "../../util/Conf.h"
#include "../../util/SVDUtil.h"
#include "../../util/Monitor.h"
#include "../../util/Calculator.h"
#include "../../structs/RedSVDMatrixRow.h"

// svd incr first done over transformed vector
class TransformSVDIncrPrune {
private:
    Monitor tt;
    double offlineTime;

    ExtendMatrix<RedSVDMatrixRow> *preprocessedP;
    Matrix *p;
    Matrix *q;
    Matrix *u;
    int k;
    int checkDim;
    int numOfNewDim;
    int numOfDim;
    vector<double> addend;
    double minValue;

    // log
#ifdef TIME_IT
    uint64_t counter1 = 0;
    uint64_t counter2 = 0;
#endif

    void transferQ(double &qNorm, double &newSVDQNorm, double &subQNorm,
                   double &sumOfCoordinate, const double *qPtr, double *newQ);
    void refine(vector<VectorElement> &heap,
                double *newQ, const double qNorm,
                const double newSVDQNorm, const double &subQNorm,
                const double sumOfQCoordinate);
public:
    TransformSVDIncrPrune(const int k, const double SIGMA, Matrix *q, Matrix *p);
    ~TransformSVDIncrPrune();
    void topK();

};

inline void TransformSVDIncrPrune::transferQ(double &qNorm, double &newSVDQNorm, double &subQNorm,
                                              double &sumOfCoordinate, const double *qPtr, double *newQ) {

    Calculator::calSingleNorm(qPtr, numOfDim, qNorm);

    newSVDQNorm = 0;
    for (int rowIndex = u->rowNum - 1; rowIndex >= 0; rowIndex--) {
        int offsetRowIndex = rowIndex + 2;
        newQ[offsetRowIndex] = 0;
        const double *uPtr = u->getRowPtr(rowIndex);

        for (int colIndex = 0; colIndex < u->colNum; colIndex++) {
            newQ[offsetRowIndex] += uPtr[colIndex] * qPtr[colIndex];
        }
        newSVDQNorm += newQ[offsetRowIndex] * newQ[offsetRowIndex];
    }
    newSVDQNorm = sqrt(newSVDQNorm);

    subQNorm = 0;
    sumOfCoordinate = 0;

    for (int colIndex = numOfNewDim-1; colIndex >= checkDim; colIndex--) {
        sumOfCoordinate += addend[colIndex - 2] * newQ[colIndex] / newSVDQNorm;
        newQ[colIndex] = (newQ[colIndex] / newSVDQNorm + addend[colIndex-2]) * 2;
        subQNorm += newQ[colIndex] * newQ[colIndex];
    }
    subQNorm = sqrt(subQNorm);

    for (int colIndex = checkDim-1; colIndex >= 2; colIndex--) {
        sumOfCoordinate += addend[colIndex - 2] * newQ[colIndex] / newSVDQNorm;
        newQ[colIndex] = (newQ[colIndex] / newSVDQNorm + addend[colIndex-2]) * 2;
    }
    newQ[0] = -1;
    newQ[1] = 0;
}

inline void TransformSVDIncrPrune::refine(vector<VectorElement> &heap,
                                           double *newQ, const double qNorm,
                                           const double newSVDQNorm, const double &subQNorm,
                                           const double sumOfQCoordinate) {
    
    for (int rowIndex = 0; rowIndex < k; rowIndex++) {
        const RedSVDMatrixRow *pRowPtr = preprocessedP->getRowPtr(rowIndex);
        double innerProduct = Calculator::innerProduct(newQ, pRowPtr->rawData, numOfNewDim);
        heap[rowIndex] = VectorElement(pRowPtr->gRowID, innerProduct);
    }

    make_heap(heap.begin(), heap.end(), greater<VectorElement>());
    double lowerBoundInNewSpace = heap.front().data;

    int gRowIDForBoundItem = heap.front().id;
    const RedSVDMatrixRow *pRowPtr = preprocessedP->getRowPtr(gRowIDForBoundItem);

    double originalLowerBound = (lowerBoundInNewSpace - pRowPtr->partialSumOfCoordinate) * 0.5 - pRowPtr->sumOfCoordinate - sumOfQCoordinate;
    originalLowerBound *= newSVDQNorm;
    originalLowerBound = isnan(originalLowerBound)?0:originalLowerBound;

    for (int rowIndex = k; rowIndex < preprocessedP->rowNum; rowIndex++) {

        const RedSVDMatrixRow *pRowPtr = preprocessedP->getRowPtr(rowIndex);

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

            if (innerProduct + pRowPtr->subVNorm * subQNorm <= lowerBoundInNewSpace) {
                continue;
            } else {

                for (int colIndex = checkDim; colIndex < numOfNewDim; colIndex++) {
                    innerProduct += newQ[colIndex] * pPtr[colIndex];
                }

#ifdef TIME_IT
                counter2++;
#endif

                if (innerProduct > lowerBoundInNewSpace) {
                    pop_heap(heap.begin(), heap.end(), greater<VectorElement>());
                    heap.pop_back();
                    heap.push_back(VectorElement(pRowPtr->gRowID, innerProduct));
                    push_heap(heap.begin(), heap.end(), greater<VectorElement>());

                    lowerBoundInNewSpace = heap.front().data;

                    int gRowIDForBoundItem = heap.front().id;
                    const RedSVDMatrixRow *pRowPtr = preprocessedP->getRowPtr(gRowIDForBoundItem);
                    originalLowerBound = (lowerBoundInNewSpace - pRowPtr->partialSumOfCoordinate) * 0.5 - pRowPtr->sumOfCoordinate - sumOfQCoordinate;
                    originalLowerBound *= newSVDQNorm;
                    originalLowerBound = isnan(originalLowerBound)?0:originalLowerBound;
                }
            }

        }
    }

}


inline TransformSVDIncrPrune::TransformSVDIncrPrune(const int k, const double SIGMA, Matrix *q, Matrix *p) {

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
    this->checkDim = SVD(P_t, p->colNum, p->rowNum, *u, *v, addend, SIGMA) + 2;

    preprocessedP = new ExtendMatrix<RedSVDMatrixRow>();
    preprocessedP->initSVDTransformExtendMatrix(*p, *v, checkDim, addend, minValue);

    delete v;

    tt.stop();
    offlineTime = tt.getElapsedTime();
}

inline TransformSVDIncrPrune::~TransformSVDIncrPrune() {

    if (u) {
        delete u;
    }

    if (preprocessedP) {
        delete preprocessedP;
    }

}

inline void TransformSVDIncrPrune::topK() {

    tt.start();

    vector<vector<VectorElement> > results(q->rowNum, vector<VectorElement>());

    double *newQ = new double[numOfNewDim];
    double subQNorm = 0;
    double newSVDQNorm = 0;
    double qNorm = 0;
    double sumOfQCoordinate = 0;

    vector<string> logStr;
    for (int qID = 0; qID < q->rowNum; qID++) {

        vector<VectorElement> &heap = results[qID];
        heap.resize(k);

        const double *qPtr = q->getRowPtr(qID);

        transferQ(qNorm, newSVDQNorm, subQNorm, sumOfQCoordinate, qPtr, newQ);

        refine(heap, newQ, qNorm, newSVDQNorm, subQNorm, sumOfQCoordinate);
    }

    for(auto s:logStr){
        cout << s << endl;
    }

    delete[] newQ;
    tt.stop();

#ifdef TIME_IT
    Logger::Log("Avg Num of p which can pass first Cauchy Schwarz inequality check: " + to_string((counter1 + 0.0)/ q
    ->rowNum));
    Logger::Log("Avg Num of p need to be calculated exactly: " + to_string((counter2 + 0.0) / q->rowNum));
#endif

    Logger::Log("preprocess time: " + to_string(offlineTime) + " secs");
    Logger::Log("online time: " + to_string(tt.getElapsedTime()) + " secs");

    if (Conf::outputResult) {
        string resultFileName = Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(k) + ".txt";
        vector<vector<VectorElement> > originalResults(q->rowNum, vector<VectorElement>(k));

        for (int qIndex = 0; qIndex < results.size(); qIndex++) {
            for (int i = 0; i < k; i++) {
                originalResults[qIndex][i] = VectorElement(results[qIndex][i].id, Calculator::innerProduct(q->getRowPtr(qIndex), p->getRowPtr(results[qIndex][i].id), numOfDim));
            }
        }
        FileUtil::outputResult(k, originalResults, resultFileName);
    }
}

#endif //TRANSFORMSVDINCRPRUNE_H
