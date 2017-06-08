#ifndef SVDINCRPRUNE_H
#define SVDINCRPRUNE_H

#include "../../util/Base.h"
#include "../../structs/Matrix.h"
#include "../../util/Conf.h"
#include "../../util/SVDUtil.h"
#include "../../util/Monitor.h"
#include "../../util/Calculator.h"
#include "../../structs/SVDMatrixRow.h"
#include "../../structs/ExtendMatrix.h"


// FEIPR-S
// only check one dimension for incremental prune
class SVDIncrPrune {

private:
    Monitor tt;
    double offlineTime;

    ExtendMatrix<SVDMatrixRow> *preprocessedP;
    Matrix *q;
    Matrix *u;
    int k;
    int checkDim;

    // log
#ifdef TIME_IT
    uint64_t counter1 = 0;
    uint64_t counter2 = 0;
    uint64_t counter3 = 0;
    uint64_t counter4 = 0;
#endif

    void transferQ(double &newQNorm, double &subQNorm, const double *qPtr, double *newQ);
    void refine(vector<VectorElement> &heap, const double *qPtr, double *newQ, const double newQNorm,
                       const double &subQNorm);
public:
    SVDIncrPrune(const int k, const double SIGMA, Matrix *q, Matrix *p);
    ~SVDIncrPrune();
    void topK();
};


inline void SVDIncrPrune::transferQ(double &newQNorm, double &subQNorm, const double *qPtr, double *newQ) {

    newQNorm = 0;
    for (int rowIndex = u->rowNum - 1; rowIndex >= 0; rowIndex--) {

        newQ[rowIndex] = 0;
        const double *uPtr = u->getRowPtr(rowIndex);

        for (int colIndex = 0; colIndex < u->colNum; colIndex++) {
            newQ[rowIndex] += uPtr[colIndex] * qPtr[colIndex];
        }

        newQNorm += newQ[rowIndex] * newQ[rowIndex];

        if (rowIndex == checkDim) {
            subQNorm = sqrt(newQNorm);
        }

    }

    newQNorm = sqrt(newQNorm);
}

inline void SVDIncrPrune::refine(vector<VectorElement> &heap, const double *qPtr, double *newQ, const double newQNorm,
                                 const double &subQNorm) {
    for (int rowIndex = 0; rowIndex < k; rowIndex++) {
        const SVDMatrixRow *pRowPtr = preprocessedP->getRowPtr(rowIndex);
        double innerProduct = Calculator::innerProduct(newQ, pRowPtr->vRawData, preprocessedP->colNum);
        heap[rowIndex] = VectorElement(pRowPtr->gRowID, innerProduct);
    }

    make_heap(heap.begin(), heap.end(), greater<VectorElement>());
    double lowerBound = heap.front().data;

    double qNorm = 0;
    Calculator::calSingleNorm(qPtr, q->colNum, qNorm);

    for (int rowIndex = k; rowIndex < preprocessedP->rowNum; rowIndex++) {

        const SVDMatrixRow *pRowPtr = preprocessedP->getRowPtr(rowIndex);
        if (pRowPtr->norm * qNorm <= lowerBound) {
            break;
        }
        else {
#ifdef TIME_IT
            counter1++;
#endif

            if (pRowPtr->vNorm * newQNorm <= lowerBound) {
                continue;
            }

#ifdef TIME_IT
            counter2++;
#endif

            const double *vPtr = pRowPtr->vRawData;

            double innerProduct = 0;
            for (int colIndex = 0; colIndex < checkDim; colIndex++) {
                innerProduct += newQ[colIndex] * vPtr[colIndex];
#ifdef TIME_IT
                counter3++;
#endif
            }

            if (innerProduct + pRowPtr->subVNorm * subQNorm <=
                lowerBound) {
                continue;
            } else {
                for (int colIndex = checkDim; colIndex < preprocessedP->colNum; colIndex++) {
                    innerProduct += newQ[colIndex] * vPtr[colIndex];
#ifdef TIME_IT
                    counter3++;
#endif
                }

#ifdef TIME_IT
                counter4++;
#endif

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
}


inline SVDIncrPrune::SVDIncrPrune(const int k, const double SIGMA, Matrix *q, Matrix *p) {

    tt.start();

    this->q = q;

    mat P_t(p->rawData, p->colNum, p->rowNum, false, true);

    vector<VectorElement> pNorms(p->rowNum);
    Calculator::calNorms(*p, pNorms);
    sort(pNorms.begin(), pNorms.end(), greater<VectorElement>());

    this->u = new Matrix();
    this->k = k;

    Matrix *v = new Matrix();
    this->checkDim = SVD(P_t, p->colNum, p->rowNum, *u, *v, SIGMA);

    vector<double> vNorms(v->rowNum);
    vector<double> subVNorms(v->rowNum);

    for (int vID = 0; vID < v->rowNum; vID++) {
        double norm = 0;
        const double *vPtr = v->getRowPtr(vID);

        for (int colIndex = v->colNum - 1; colIndex >= 0; colIndex--) {
            norm += vPtr[colIndex] * vPtr[colIndex];
            if (colIndex == checkDim) {
                subVNorms[vID] = sqrt(norm);
            }
        }
        norm = sqrt(norm);
        vNorms[vID] = norm;
    }

    preprocessedP = new ExtendMatrix<SVDMatrixRow>();
    preprocessedP->initSVDExtendMatrix(*p, *v, pNorms, vNorms, subVNorms);

    delete v;

    tt.stop();
    offlineTime = tt.getElapsedTime();
}

inline SVDIncrPrune::~SVDIncrPrune() {

    if (u) {
        delete u;
    }

    if (preprocessedP) {
        delete preprocessedP;
    }

}

inline void SVDIncrPrune::topK() {

    tt.start();

    vector<vector<VectorElement> > results(q->rowNum, vector<VectorElement>());

    double *newQ = new double[q->colNum];
    double subQNorm = 0;
    double newQNorm = 0;

    for (int qID = 0; qID < q->rowNum; qID++) {

        // step 1: transfer q

        vector<VectorElement> &heap = results[qID];
        heap.resize(k);

        const double *qPtr = q->getRowPtr(qID);
        transferQ(newQNorm, subQNorm, qPtr, newQ);

        // step 3: refine
        refine(heap, qPtr, newQ, newQNorm, subQNorm);

    }

    delete[] newQ;
    tt.stop();

#ifdef TIME_IT
    Logger::Log("Avg Num of p which can pass first Cauchy Schwarz inequality check: " + to_string((counter1 + 0.0)/ q
    ->rowNum));
    Logger::Log("Avg Num of p which can pass second Cauchy Schwarz inequality check: " + to_string((counter2 + 0.0) / q->rowNum));
    Logger::Log("Avg Num of Dimension need to be calculated: " + to_string((counter3 + 0.0) / counter2));
    Logger::Log("Avg Num of p need to be calculated exactly: " + to_string((counter4 + 0.0) / q->rowNum));
#endif

    Logger::Log("preprocess time: " + to_string(offlineTime) + " secs");
    Logger::Log("online time: " + to_string(tt.getElapsedTime()) + " secs");

    if (Conf::outputResult) {
        string resultFileName =
                Conf::resultPathPrefix + "-" + Conf::dataset +  "-" + Conf::algName + "-" + to_string(Conf::k) + "-" +
                to_string(Conf::SIGMA) + ".txt";
        FileUtil::outputResult(k, results, resultFileName);
    }
}

#endif //SVDINCRPRUNE_H
