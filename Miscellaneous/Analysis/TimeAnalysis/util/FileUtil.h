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

    void outputFlag(int *flag, int size) {
        string fileName = "flag.txt";
        ofstream file(fileName.c_str());

        for (int qIndex = 0; qIndex < size; qIndex++) {
            file << flag[qIndex] << endl;
        }
        file.close();
    }

    void outputResult(vector<vector<VectorElement> > &results, string &resultPath) {
        createFolder(resultPath);
        ofstream file(resultPath.c_str());

        for (int qIndex = 0; qIndex < results.size(); qIndex++) {
            for (int kIndex = 0; kIndex < results[qIndex].size(); kIndex++) {
                file << qIndex << " " << results[qIndex][kIndex].id << " " <<
                        results[qIndex][kIndex].data << endl;
            }
        }
        file.close();
    }

    void readFlag(int *flags) {

        string line;
        string dataFilePath = "flag.txt";

        ifstream fin(dataFilePath.c_str());

        int index = 0;
        while (getline(fin, line)) {
            if (line.length() == 0) {
                continue;
            }
            flags[index] = atoi(line.c_str());
            index++;
        }
    }

}
#endif //FILEUTIL_H
