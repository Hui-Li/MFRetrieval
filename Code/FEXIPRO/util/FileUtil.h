#ifndef FILEUTIL_H
#define FILEUTIL_H

#include "Base.h"
#include "../structs/VectorElement.h"

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

    void outputResult(const int k, vector<vector<VectorElement> > &results, string &resultPath) {
        createFolder(resultPath);
        ofstream file(resultPath.c_str());

        if(k>1){
            vector<VectorElement> toSort(k);

            for (int qIndex = 0; qIndex < results.size(); qIndex++) {
                for (int kIndex = 0; kIndex < results[qIndex].size(); kIndex++) {
                    toSort[kIndex] = results[qIndex][kIndex];
                }

                sort(toSort.begin(), toSort.end(), greater<VectorElement>());
                for (int kIndex = 0; kIndex < k; kIndex++) {
                    file << qIndex << " " << toSort[kIndex].id << " " << toSort[kIndex].data << endl;
                }
            }
        } else {
            for (int qIndex = 0; qIndex < results.size(); qIndex++) {
                for (int kIndex = 0; kIndex < results[qIndex].size(); kIndex++) {
                    file << qIndex << " " << results[qIndex][kIndex].id << " " << results[qIndex][kIndex].data << endl;
                }
            }
        }

        file.close();
    }

    void outputResult(const int qNum, const int k, VectorElement *results, string &resultPath) {
        createFolder(resultPath);
        ofstream file(resultPath.c_str());

        if(k>1){
            vector<VectorElement> toSort(k);
            for (int qIndex = 0; qIndex < qNum; qIndex++) {
                for (int kIndex = 0; kIndex < k; kIndex++) {
                    toSort[kIndex] = results[qIndex*k+kIndex];
                }

                sort(toSort.begin(), toSort.end(), greater<VectorElement>());
                for (int kIndex = 0; kIndex < k; kIndex++) {
                    file << qIndex << " " << toSort[kIndex].id << " " << toSort[kIndex].data << endl;
                }
            }
        }else{
            for (int qIndex = 0; qIndex < qNum; qIndex++) {
                for (int kIndex = 0; kIndex < k; kIndex++) {
                    file << qIndex << " " << results[qIndex*k+kIndex].id << " " << results[qIndex*k+kIndex].data << endl;
                }
            }
        }

        file.close();
    }
}
#endif //FILEUTIL_H
