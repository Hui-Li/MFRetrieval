#ifndef TRANSFORMINCRPRUNE_H
#define TRANSFORMINCRPRUNE_H

#include "../../util/Base.h"
#include "../../util/TransformUtil.h"
#include "../../util/Monitor.h"
#include "../../structs/Matrix.h"
#include "../../util/Calculator.h"
#include "../../structs/ExtendMatrixRow.h"
#include "../../structs/ExtendMatrix.h"

class TransformIncrPrune {

private:
    Monitor tt;
    double offlineTime;
    Matrix *q;
    Matrix *p;
    ExtendMatrix<ExtendMatrixRow> *preprocessedP;
    vector<double> addend;
    int checkDim;
    int k;
    int numOfNewDim;
    int numOfDim;

public:
    TransformIncrPrune(const int k, const int checkDim, Matrix *q, Matrix *p);
    ~TransformIncrPrune();
    void transferQ(double &newQNorm, double &subQNorm, const double *qPtr, double *newQ);
    void topK();
};

inline TransformIncrPrune::TransformIncrPrune(const int k, const int checkDim, Matrix *q, Matrix *p){

    tt.start();
    this->k = k;
    this->checkDim = checkDim;
    this->q = q;
    this->p = p;
    this->numOfDim = q->colNum;
    this->numOfNewDim = q->colNum + 2;

    vector<VectorElement> pNorms(p->rowNum);
    double maxPNorm;
    Calculator::calNorms(*p, pNorms, maxPNorm);

    // ToDo: why different addend have the same result?
    addend.resize(q->colNum);
    for (int i = 0; i < q->colNum; i++) {
        // decay function
        addend[i] = 10000 * pow(0.8, i);
        addend[i] = std::max(addend[i], 1.0);
        addend[i] = q->colNum - i;
//        addend[i] = 1;
    }

    Matrix transformedP;
    vector<VectorElement> sortedNewPNorm;
    vector<double> unsortedSubNorm;

    TransformUtil::makePPositive(*p, transformedP, maxPNorm, checkDim, pNorms, sortedNewPNorm, unsortedSubNorm, addend);

    preprocessedP = new ExtendMatrix<ExtendMatrixRow>();
    preprocessedP->initExtendMatrix(transformedP, checkDim, sortedNewPNorm, unsortedSubNorm);

    tt.stop();
    offlineTime = tt.getElapsedTime();
}

inline TransformIncrPrune::~TransformIncrPrune(){
    if(preprocessedP) {
        delete preprocessedP;
    }
}

inline void TransformIncrPrune::transferQ(double &newQNorm, double &subQNorm, const double *qPtr, double *newQ) {

    newQ[0] = -1;
    newQ[1] = 2;

    newQNorm = 0;
    double qNorm;
    Calculator::calSingleNorm(qPtr, numOfDim, qNorm);

    for (int colIndex = numOfNewDim-1; colIndex >= 2; colIndex--) {

        newQ[colIndex] = (qPtr[colIndex-2] / qNorm + addend[colIndex-2]) * 2;
        newQNorm += newQ[colIndex] * newQ[colIndex];
        if (colIndex == checkDim) {
            subQNorm = sqrt(newQNorm);
        }

    }

    newQNorm = sqrt(newQNorm + 5) - 1;
}

inline void TransformIncrPrune::topK(){

#ifdef TIME_IT
    uint64_t counter1 = 0;
    uint64_t counter2 = 0;
#endif

    tt.start();

    vector<vector<VectorElement> > results(q->rowNum, vector<VectorElement>());

    double newQNorm;
    double subQNorm;
    double *newQ = new double[q->colNum + 2];

    for (int qID = 0; qID < q->rowNum; qID++) {

        const double *qRow = q->getRowPtr(qID);
        transferQ(newQNorm, subQNorm, qRow, newQ);

        vector<VectorElement> &heap = results[qID];
        heap.resize(k);

        for (int pIndex = 0; pIndex < k; pIndex++) {
            const ExtendMatrixRow * pRow = preprocessedP->getRowPtr(pIndex);
            double innerProduct = Calculator::innerProduct(newQ, pRow->rawData, numOfNewDim);
            heap[pIndex] = VectorElement(pRow->gRowID, innerProduct);
        }

        make_heap(heap.begin(), heap.end(), greater<VectorElement>());
        double lowerBound = heap.front().data;

        for (int pIndex = k; pIndex < preprocessedP->rowNum; pIndex++) {

            const ExtendMatrixRow * pRow = preprocessedP->getRowPtr(pIndex);

            if (pRow->norm * newQNorm <= lowerBound) {
                break;
            }
            else {

#ifdef TIME_IT
                counter1++;
#endif
                double innerProduct = 0;
                const double *pPtr = pRow->rawData;
                for (int colIndex = 0; colIndex < checkDim; colIndex++) {
                    innerProduct += newQ[colIndex] * pPtr[colIndex];
                }

                if(innerProduct + subQNorm * pRow->subNorm < lowerBound){
                    continue;
                }

#ifdef TIME_IT
                counter2++;
#endif

                for (int colIndex = checkDim; colIndex < numOfNewDim; colIndex++) {
                    innerProduct += newQ[colIndex] * pPtr[colIndex];
                }

                if (innerProduct > lowerBound) {
                    pop_heap(heap.begin(), heap.end(), greater<VectorElement>());
                    heap.pop_back();
                    heap.push_back(VectorElement(pRow->gRowID, innerProduct));
                    push_heap(heap.begin(), heap.end(), greater<VectorElement>());
                    lowerBound = heap.front().data;
                }
            }
        }
//        cout << heap[0].id << "," << heap[0].data << "," << Calculator::innerProduct(q->getRowPtr(qID), p->getRowPtr(heap[0].id), numOfDim) << endl;
//        exit(0);
    }

    tt.stop();

    Logger::Log("preprocess time: " + to_string(offlineTime) + " secs");
    Logger::Log("online time: " + to_string(tt.getElapsedTime()) + " secs");

#ifdef TIME_IT
    Logger::Log("Avg Num of p which can pass first Cauchy Schwarz inequality check: " + to_string((counter1 + 0.0)/ q
        ->rowNum));
    Logger::Log("Avg Num of p which can pass incremental prune: " + to_string((counter2 + 0.0) / q->rowNum));
#endif

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
#endif //TRANSFORMINCRPRUNE_H
