#ifndef ARMANAIVE_H
#define ARMANAIVE_H

#include "../util/Base.h"
#include "../util/FileUtil.h"
#include "../util/Conf.h"
#include "../util/Monitor.h"
#include "../util/Logger.h"
#include "../structs/VectorElement.h"

namespace ArmaNaive {

    vector<vector<VectorElement> > results;
    int index = 0;
    int size;

    void naive() {

        int qRowNum;
        int qColNum;

        FileUtil::getSize(Conf::qDataPath, qRowNum, qColNum);
        double *qRawData = new double[qRowNum * qColNum];
        FileUtil::readData(Conf::qDataPath, qRawData, qColNum);

        mat Q(qRawData, qColNum, qRowNum, false, true);
        Q = Q.t();

//        cout << Q << endl;

        int pRowNum;
        int pColNum;

        FileUtil::getSize(Conf::pDataPath, pRowNum, pColNum);
        double *pRawData = new double[pRowNum * pColNum];
        FileUtil::readData(Conf::pDataPath, pRawData, pColNum);

        mat P_t(pRawData, pColNum, pRowNum, false, true);

//        cout << P_t << endl;

        Logger::Log("q: " + to_string(qRowNum) + "," + to_string(qColNum));
        Logger::Log("p: " + to_string(pRowNum) + "," + to_string(pColNum));

        Monitor t;
        t.start();

        mat dot_product = Q * P_t;

//        cout << dot_product << endl;

        size = std::min(Conf::k, qColNum);
        results.resize(qRowNum, vector<VectorElement>(size));

        dot_product.each_row([](const rowvec &a) {

            uvec sortedIds = sort_index(a, "descend");
            colvec topk = a.elem(sortedIds);
            vector<VectorElement> &result = results[index];

            for (int i = 0; i < size; i++) {
                result[i].id = sortedIds(i);
                result[i].data = topk(i);
            }

            index++;

        });

//        for (int rowIndex = 0; rowIndex < results.size(); rowIndex++) {
//            vector<VectorElement> &result = results[rowIndex];
//            for (int colIndex = 0; colIndex < result.size(); colIndex++) {
//                cout << rowIndex << "," << result[colIndex].id << "," << result[colIndex].data << endl;
//            }
//        }

        t.stop();

        delete[] qRawData;
        delete[] pRawData;

        Logger::Log("time: " + to_string(t.getElapsedTime()) + " secs");

        if (Conf::outputResult) {
            string resultFileName =
                    Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) +
                    ".txt";
            FileUtil::outputResult(results, resultFileName);
        }

    }

    void batchNaive() {

        int qRowNum;
        int qColNum;

        FileUtil::getSize(Conf::qDataPath, qRowNum, qColNum);
        double *qRawData = new double[qRowNum * qColNum];
        FileUtil::readData(Conf::qDataPath, qRawData, qColNum);

        mat Q(qRawData, qColNum, qRowNum, false, true);
        Q = Q.t();

//        cout << Q << endl;

        int pRowNum;
        int pColNum;

        FileUtil::getSize(Conf::pDataPath, pRowNum, pColNum);
        double *pRawData = new double[pRowNum * pColNum];
        FileUtil::readData(Conf::pDataPath, pRawData, pColNum);

        mat P_t(pRawData, pColNum, pRowNum, false, true);

//        cout << P_t << endl;

        Logger::Log("q: " + to_string(qRowNum) + "," + to_string(qColNum));
        Logger::Log("p: " + to_string(pRowNum) + "," + to_string(pColNum));

//        cout << Q * P_t << endl;

        Monitor t;
        t.start();

        int batch = qRowNum / Conf::batchSize;

        size = std::min(Conf::k, qColNum);
        results.resize(qRowNum, vector<VectorElement>(size));

        for (int batchID = 0; batchID < batch; batchID++) {
            mat dot_product = Q.rows(batchID * Conf::batchSize, (batchID + 1) * Conf::batchSize - 1) * P_t;

            dot_product.each_row([](const rowvec &a) {

                uvec sortedIds = sort_index(a, "descend");
                colvec topk = a.elem(sortedIds);
                vector<VectorElement> &result = results[index];

                for (int i = 0; i < size; i++) {
                    result[i].id = sortedIds(i);
                    result[i].data = topk(i);
                }

                index++;

            });
        }

        if (qRowNum % Conf::batchSize != 0) {

            mat dot_product = Q.rows(batch * Conf::batchSize, qRowNum) * P_t;

            dot_product.each_row([](const rowvec &a) {

                uvec sortedIds = sort_index(a, "descend");
                colvec topk = a.elem(sortedIds);
                vector<VectorElement> &result = results[index];

                for (int i = 0; i < size; i++) {
                    result[i].id = sortedIds(i);
                    result[i].data = topk(i);
                }

                index++;

            });

            Logger::Log("batchNum: " + to_string(batch + 1));

        } else {
            Logger::Log("batchNum: " + to_string(batch));
        }

//        for (int rowIndex = 0; rowIndex < results.size(); rowIndex++) {
//            vector<VectorElement> &result = results[rowIndex];
//            for (int colIndex = 0; colIndex < result.size(); colIndex++) {
//                cout << rowIndex << "," << result[colIndex].id << "," << result[colIndex].data << endl;
//            }
//        }

        t.stop();

        delete[] qRawData;
        delete[] pRawData;

        Logger::Log("time: " + to_string(t.getElapsedTime()) + " secs");

        if (Conf::outputResult) {
            string resultFileName =
                    Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) +
                    ".txt";
            FileUtil::outputResult(results, resultFileName);
        }
    }
}

#endif //ARMANAIVE_H
