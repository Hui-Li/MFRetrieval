#ifndef FILEUTIL_H
#define FILEUTIL_H

#include "Base.h"
#include "../structs/VectorElement.h"
#include "../structs/Triplet.h"
#include "../structs/Matrix.h"

namespace FileUtil {

    void createFolder(string folderPath) {
        for (int i = 0; i < folderPath.length(); ++i) {
            if (folderPath[i] == '\\') {
                folderPath[i] = '/';
            }
            if (folderPath[i] == '/') {
                string folder = folderPath.substr(0, i);
#ifdef _WIN32
                _mkdir(folder.c_str());
#else
                mkdir(folder.c_str(), 0777);
#endif
            }
        }
    }

    void outputResult(const int k, Matrix &P, Matrix &Q, VectorElement *results, string &resultPath) {
        createFolder(resultPath);
        ofstream file(resultPath.c_str());

        if (k > 1) {
            vector<Triplet> toSort(k);
            for (int qIndex = 0; qIndex < Q.rowNum; qIndex++) {
                for (int kIndex = 0; kIndex < k; kIndex++) {
                    int pIndex = results[qIndex * k + kIndex].id;
                    toSort[kIndex].id = pIndex;
                    toSort[kIndex].value1 = inner_product(Q[qIndex], P[pIndex], Q.colNum);
                    toSort[kIndex].value2 = results[qIndex * k + kIndex].data;
                }

                sort(toSort.begin(), toSort.end(), std::greater<Triplet>());
                for (int kIndex = 0; kIndex < k; kIndex++) {
                    file << qIndex << "," << toSort[kIndex].id << "," << toSort[kIndex].value1 << "," << toSort[kIndex].value2 << endl;
                }
            }
        } else {
            for (int qIndex = 0; qIndex < Q.rowNum; qIndex++) {
                VectorElement* sub_result = results + qIndex;
                double real_product = inner_product(Q[qIndex], P[sub_result->id], Q.colNum);
                file << qIndex << "," << sub_result->id << "," << real_product
                     << "," << sub_result->data << endl;
            }
        }

        file.close();
    }

    void outputNaiveResult(const int k, Matrix &P, Matrix &Q, VectorElement *results, string &resultPath) {
        createFolder(resultPath);
        ofstream file(resultPath.c_str());

        if (k > 1) {
            for (int qIndex = 0; qIndex < Q.rowNum; qIndex++) {
                VectorElement *sub_results = results + qIndex * k;
                sort(sub_results, sub_results + k, std::greater<VectorElement>());
                for (int kIndex = 0; kIndex < k; kIndex++) {
                    file << qIndex << "," << sub_results[kIndex].id << "," << sub_results[kIndex].data << endl;
                }
            }
        } else {
            for (int qIndex = 0; qIndex < Q.rowNum; qIndex++) {
                VectorElement *sub_result = results + qIndex;
                file << qIndex << "," << sub_result->id << "," << sub_result->data << endl;
            }
        }

        file.close();
    }
}
#endif //FILEUTIL_H
