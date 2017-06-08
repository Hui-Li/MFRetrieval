#ifndef MKLNAIVE_H
#define MKLNAIVE_H

#include <stdio.h>
#include "../util/Base.h"
#include "mkl.h"

namespace MKLNaive {

    // https://software.intel.com/en-us/node/529735
    void naive() {

        if (Conf::singleThread) {
            mkl_set_num_threads(1);
        }

        // step 1: get all dot products
        int qRowNum;
        int qColNum;

        FileUtil::getSize(Conf::qDataPath, qRowNum, qColNum);
        double *qRawData = (double *) mkl_malloc(qRowNum * qColNum * sizeof(double), 64);
        FileUtil::readData(Conf::qDataPath, qRawData, qColNum);

        int pRowNum;
        int pColNum;

        FileUtil::getSize(Conf::pDataPath, pRowNum, pColNum);
        double *pRawData = (double *) mkl_malloc(pRowNum * pColNum * sizeof(double), 64);
        FileUtil::readDataT(Conf::pDataPath, pRawData, pColNum);

        double alpha = 1.0, beta = 0.0;

        Logger::Log("q: " + to_string(qRowNum) + "," + to_string(qColNum));
        Logger::Log("p: " + to_string(pRowNum) + "," + to_string(pColNum));

        vector<vector<VectorElement> > results(qRowNum);

        Monitor t;
        t.start();

        double *dotProduct = (double *) mkl_malloc(qRowNum * pRowNum * sizeof(double), 64);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    qRowNum, pRowNum, qColNum, alpha, qRawData, qColNum, pRawData, pRowNum, beta, dotProduct, pRowNum);

//        for (int i = 0; i < qRowNum; i++) {
//            for (int j = 0; j < pRowNum; j++) {
//                cout << dotProduct[i * pRowNum + j] << " ";
//            }
//            cout << endl;
//        }

        mkl_free(qRawData);
        mkl_free(pRawData);

        // step 2: get top-k
        int tmpRow;
        for (int qID = 0; qID < qRowNum; qID++) {
            vector<VectorElement> &heap = results[qID];
            heap.resize(Conf::k);

            tmpRow = qID * pRowNum;
            for (int pID = 0; pID < Conf::k; pID++) {
                heap[pID] = VectorElement(pID, dotProduct[tmpRow + pID]);
            }

            make_heap(heap.begin(), heap.end(), std::greater<VectorElement>());
            double lowerBound = heap.front().data;

            for (int pID = Conf::k; pID < pRowNum; pID++) {
                double innerProduct = dotProduct[tmpRow + pID];

                if (innerProduct > lowerBound) {
                    pop_heap(heap.begin(), heap.end(), std::greater<VectorElement>());
                    heap.pop_back();
                    heap.push_back(VectorElement(pID, innerProduct));
                    push_heap(heap.begin(), heap.end(), std::greater<VectorElement>());
                    lowerBound = heap.front().data;
                }
            }
        }

        t.stop();

        Logger::Log("time: " + to_string(t.getElapsedTime()) + " secs");

        if (Conf::outputResult) {
            string resultFileName =
                    Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k);

            if (Conf::singleThread) {
                resultFileName += "-singlethread.txt";
            } else {
                resultFileName += ".txt";
            }

            FileUtil::outputResult(results, resultFileName);
        }

        mkl_free(dotProduct);

//        for (int rowIndex = 0; rowIndex < results.size(); rowIndex++) {
//            vector<VectorElement> &result = results[rowIndex];
//            for (int colIndex = 0; colIndex < result.size(); colIndex++) {
//                cout << rowIndex << "," << result[colIndex].id << "," << result[colIndex].data << endl;
//            }
//        }

    }

    void batchNaive() {

        if (Conf::singleThread) {
            mkl_set_num_threads(1);
        }

        int qRowNum;
        int qColNum;

        FileUtil::getSize(Conf::qDataPath, qRowNum, qColNum);
        double *qRawData = (double *) mkl_malloc(qRowNum * qColNum * sizeof(double), 64);
        FileUtil::readData(Conf::qDataPath, qRawData, qColNum);

        int pRowNum;
        int pColNum;

        FileUtil::getSize(Conf::pDataPath, pRowNum, pColNum);
        double *pRawData = (double *) mkl_malloc(pRowNum * pColNum * sizeof(double), 64);
        FileUtil::readDataT(Conf::pDataPath, pRawData, pColNum);

        double alpha = 1.0, beta = 0.0;

        Logger::Log("q: " + to_string(qRowNum) + "," + to_string(qColNum));
        Logger::Log("p: " + to_string(pRowNum) + "," + to_string(pColNum));

        vector<vector<VectorElement> > results(qRowNum);

        Monitor t;
        t.start();

        int batch = qRowNum / Conf::batchSize;
        int skipSizeForOneBatch = Conf::batchSize * qColNum;

        double *dotProduct = (double *) mkl_malloc(Conf::batchSize * pRowNum * sizeof(double), 64);
        int tmpRow;

        for (int batchID = 0; batchID < batch; batchID++) {

            double *subQRawData = &qRawData[batchID * skipSizeForOneBatch];

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        Conf::batchSize, pRowNum, qColNum, alpha, subQRawData, qColNum, pRawData, pRowNum, beta, dotProduct, pRowNum);

            for (int localQID = 0; localQID < Conf::batchSize; localQID++) {
                int globalQID = batchID * Conf::batchSize + localQID;
                vector<VectorElement> &heap = results[globalQID];
                heap.resize(Conf::k);

                tmpRow = localQID * pRowNum;
                for (int pID = 0; pID < Conf::k; pID++) {
                    heap[pID] = VectorElement(pID, dotProduct[tmpRow + pID]);
                }

                make_heap(heap.begin(), heap.end(), std::greater<VectorElement>());
                double lowerBound = heap.front().data;

                for (int pID = Conf::k; pID < pRowNum; pID++) {
                    double innerProduct = dotProduct[tmpRow + pID];

                    if (innerProduct > lowerBound) {
                        pop_heap(heap.begin(), heap.end(), std::greater<VectorElement>());
                        heap.pop_back();
                        heap.push_back(VectorElement(pID, innerProduct));
                        push_heap(heap.begin(), heap.end(), std::greater<VectorElement>());
                        lowerBound = heap.front().data;
                    }
                }
            }

        }

        int remainQNum = qRowNum % Conf::batchSize;

        if (remainQNum != 0) {

            int start = batch * Conf::batchSize;

            double *subQRawData = &qRawData[batch * skipSizeForOneBatch];

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        remainQNum, pRowNum, qColNum, alpha, subQRawData, qColNum, pRawData, pRowNum, beta, dotProduct, pRowNum);

            for (int localQID = 0; localQID < remainQNum; localQID++) {
                int globalQID = start + localQID;
                vector<VectorElement> &heap = results[globalQID];
                heap.resize(Conf::k);

                tmpRow = localQID * pRowNum;
                for (int pID = 0; pID < Conf::k; pID++) {
                    heap[pID] = VectorElement(pID, dotProduct[tmpRow + pID]);
                }

                make_heap(heap.begin(), heap.end(), std::greater<VectorElement>());
                double lowerBound = heap.front().data;

                for (int pID = Conf::k; pID < pRowNum; pID++) {
                    double innerProduct = dotProduct[tmpRow + pID];

                    if (innerProduct > lowerBound) {
                        pop_heap(heap.begin(), heap.end(), std::greater<VectorElement>());
                        heap.pop_back();
                        heap.push_back(VectorElement(pID, innerProduct));
                        push_heap(heap.begin(), heap.end(), std::greater<VectorElement>());
                        lowerBound = heap.front().data;
                    }
                }
            }

            Logger::Log("batchNum: " + to_string(batch + 1));

        } else {
            Logger::Log("batchNum: " + to_string(batch));
        }

        t.stop();

        Logger::Log("time: " + to_string(t.getElapsedTime()) + " secs");

        if (Conf::outputResult) {
            string resultFileName =
                    Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k);

            if (Conf::singleThread) {
                resultFileName += "-singlethread.txt";
            } else {
                resultFileName += ".txt";
            }

            FileUtil::outputResult(results, resultFileName);
        }

        mkl_free(qRawData);
        mkl_free(pRawData);
        mkl_free(dotProduct);

//        for (int rowIndex = 0; rowIndex < results.size(); rowIndex++) {
//            vector<VectorElement> &result = results[rowIndex];
//            for (int colIndex = 0; colIndex < result.size(); colIndex++) {
//                cout << rowIndex << "," << result[colIndex].id << "," << result[colIndex].data << endl;
//            }
//        }
    }

}

#endif //MKLNAIVE_H
