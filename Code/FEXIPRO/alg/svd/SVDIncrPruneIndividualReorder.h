#ifndef SVDINCRPRUNEINDIVIDUALREORDER_H
#define SVDINCRPRUNEINDIVIDUALREORDER_H

#include "../../util/Base.h"
#include "../../structs/Matrix.h"
#include "../../util/Conf.h"
#include "../../util/SVDUtil.h"
#include "../../util/Monitor.h"
#include "../../util/Calculator.h"
#include "../../structs/SVDMatrixRowIndividualReorder.h"
#include "../../structs/ExtendMatrix.h"


// FEIPR-S with individual reorder
// only check one dimension for incremental prune
class SVDIncrPruneIndividualReorder {

private:
    Monitor tt;
    double offlineTime;

    ExtendMatrix<SVDMatrixRowIndividualReorder> *preprocessedP;
    Matrix *q;
    Matrix *u;
    int k;
    int checkDim;
    int dim;
    double SIGMA;

    // log
#ifdef TIME_IT
    uint64_t counter1 = 0;
    uint64_t counter2 = 0;
#endif

    void transferQ(double &newQNorm, double &subQNorm, const double *qPtr, double *newQ, int &checkDim);
    void refine(vector<VectorElement> &heap, const double *qPtr, double *newQ, const double newQNorm,
                const double &subQNorm, const int checkDim);
public:
    SVDIncrPruneIndividualReorder(const int k, const double SIGMA, Matrix *q, Matrix *p);
    ~SVDIncrPruneIndividualReorder();
    void topK();
};


inline void SVDIncrPruneIndividualReorder::transferQ(double &newQNorm, double &subQNorm, const double *qPtr, double *newQ, int &checkDim) {

    newQNorm = 0;

    double *sum = new double[this->dim];

    for (int rowIndex = 0; rowIndex < u->rowNum; rowIndex++) {

        newQ[rowIndex] = 0;
        const double *uPtr = u->getRowPtr(rowIndex);

        for (int colIndex = 0; colIndex < u->colNum; colIndex++) {
            newQ[rowIndex] += uPtr[colIndex] * qPtr[colIndex];
        }

        newQNorm += newQ[rowIndex] * newQ[rowIndex];
        sum[rowIndex] = newQNorm;
    }

    for (int colIndex = 0; colIndex < dim; colIndex++) {
        if(sum[colIndex] / newQNorm >= SIGMA){
            checkDim = std::min(colIndex + 1, dim);
            subQNorm = sqrt(newQNorm - sum[colIndex]);
            break;
        }
    }

    newQNorm = sqrt(newQNorm);
    delete [] sum;
}

inline void SVDIncrPruneIndividualReorder::refine(vector<VectorElement> &heap, const double *qPtr, double *newQ, const double newQNorm,
                                 const double &subQNorm, const int checkDim) {
    for (int rowIndex = 0; rowIndex < k; rowIndex++) {
        const SVDMatrixRowIndividualReorder *pRowPtr = preprocessedP->getRowPtr(rowIndex);
        double innerProduct = Calculator::innerProduct(newQ, pRowPtr->vRawData, preprocessedP->colNum);
        heap[rowIndex] = VectorElement(pRowPtr->gRowID, innerProduct);
    }

    make_heap(heap.begin(), heap.end(), greater<VectorElement>());
    double lowerBound = heap.front().data;

    double qNorm = 0;
    Calculator::calSingleNorm(qPtr, q->colNum, qNorm);

    for (int rowIndex = k; rowIndex < preprocessedP->rowNum; rowIndex++) {

        const SVDMatrixRowIndividualReorder *pRowPtr = preprocessedP->getRowPtr(rowIndex);
        if (pRowPtr->norm * qNorm <= lowerBound) {
            break;
        }
        else {
#ifdef TIME_IT
            counter1++;
#endif

//            if (pRowPtr->vNorm * newQNorm <= lowerBound) {
//                continue;
//            }

//#ifdef TIME_IT
//            counter2++;
//#endif

            const double *vPtr = pRowPtr->vRawData;

            double innerProduct = 0;
            for (int colIndex = 0; colIndex < checkDim; colIndex++) {
                innerProduct += newQ[colIndex] * vPtr[colIndex];
            }

            if (innerProduct + pRowPtr->subVNorm[checkDim] * subQNorm <=
                lowerBound) {
                continue;
            } else {
                for (int colIndex = checkDim; colIndex < preprocessedP->colNum; colIndex++) {
                    innerProduct += newQ[colIndex] * vPtr[colIndex];

                }

#ifdef TIME_IT
                counter2++;
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


inline SVDIncrPruneIndividualReorder::SVDIncrPruneIndividualReorder(const int k, const double SIGMA, Matrix *q, Matrix *p) {

    tt.start();

    this->q = q;
    this->SIGMA = SIGMA;

    mat P_t(p->rawData, p->colNum, p->rowNum, false, true);

    vector<VectorElement> pNorms(p->rowNum);
    Calculator::calNorms(*p, pNorms);
    sort(pNorms.begin(), pNorms.end(), greater<VectorElement>());

    this->u = new Matrix();
    this->k = k;
    this->dim = p->colNum;

    Matrix *v = new Matrix();
    double *vSubNorms = new double[p->rowNum * p->colNum];
    SVD(P_t, p->colNum, p->rowNum, *u, *v, vSubNorms);
    preprocessedP = new ExtendMatrix<SVDMatrixRowIndividualReorder>();
    preprocessedP->initSVDExtendMatrix(*p, *v, pNorms, vSubNorms);

    delete v;
    delete [] vSubNorms;
    tt.stop();
    offlineTime = tt.getElapsedTime();
}

inline SVDIncrPruneIndividualReorder::~SVDIncrPruneIndividualReorder() {

    if (u) {
        delete u;
    }

    if (preprocessedP) {
        delete preprocessedP;
    }

}

inline void SVDIncrPruneIndividualReorder::topK() {

    tt.start();

    vector<vector<VectorElement> > results(q->rowNum, vector<VectorElement>());

    double *newQ = new double[q->colNum];
    double subQNorm = 0;
    double newQNorm = 0;
    int checkDim = 0;

    for (int qID = 0; qID < q->rowNum; qID++) {

        // step 1: transfer q

        vector<VectorElement> &heap = results[qID];
        heap.resize(k);

        const double *qPtr = q->getRowPtr(qID);
        transferQ(newQNorm, subQNorm, qPtr, newQ, checkDim);

        // step 3: refine
        refine(heap, qPtr, newQ, newQNorm, subQNorm, checkDim);

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
        string resultFileName =
                Conf::resultPathPrefix + "-" + Conf::dataset +  "-" + Conf::algName + "-" + to_string(Conf::k) + "-" +
                to_string(Conf::SIGMA) + ".txt";
        FileUtil::outputResult(k, results, resultFileName);
    }
}

#endif //SVDINCRPRUNEINDIVIDUALREORDER_H
