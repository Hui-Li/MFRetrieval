#ifndef FILEUTIL_H
#define FILEUTIL_H

#include "Base.h"

namespace FileUtil{

    void read(const string dataFilePath, int &rowNum, int &colNum, double *&rawData){
        vector <string> lines;
        string line;

        ifstream fin(dataFilePath.c_str());

        rowNum = 0;
        colNum = 0;

        while (getline(fin, line)) {
            if (line.length() == 0) {
                continue;
            }
            lines.push_back(line);
        }

        if (lines.size() == 0) {
            return;
        }

        fin.close();

        stringstream test(lines[0]);
        double tempValue;
        while (test >> tempValue) {
            if (test.peek() == ',') {
                test.ignore();
                colNum++;
            }
        }
        colNum++;
        rowNum = lines.size();

        rawData = new double[rowNum * colNum];
        int colIndex = 0;

        for (int rowIndex = 0; rowIndex < lines.size(); rowIndex++) {
            stringstream ss(lines[rowIndex]);
            colIndex = 0;
            while (ss >> rawData[rowIndex * colNum + colIndex]) {
                if (ss.peek() == ',') {
                    ss.ignore();
                    colIndex++;
                }
            }
        }

    }

    void outputResult(const int k, const int d, const int qNum, const double *q, const double *p, VectorElement *results, string &resultPath) {

        ofstream file(resultPath.c_str());

        for (int qIndex = 0; qIndex < qNum; qIndex++) {
            for (int i = 0; i < k; i++) {
                const double *qPtr = q + qIndex * d;
                const double *pPtr = p + results[qIndex * k + i].id * d;
                double ip = std::inner_product(qPtr, qPtr + d, pPtr, 0.0);
                file << qIndex << "," << results[qIndex * k + i].id << "," << ip << endl;
            }
        }
        file.close();
    }

    void outputResult(const int k, const int d, const int qNum, const double *q, const double *p, VectorElement *results, VectorElement *naiveResults, string &resultPath) {

        ofstream file(resultPath.c_str());

        for (int qIndex = 0; qIndex < qNum; qIndex++) {
            for (int i = 0; i < k; i++) {
                const double *qPtr = q + qIndex * d;
                const double *pPtr = p + results[qIndex * k + i].id * d;
                double ip = std::inner_product(qPtr, qPtr + d, pPtr, 0.0);
                file << qIndex << "," << results[qIndex * k + i].id << "," << ip << "," << naiveResults[qIndex * k + i].id << "," << naiveResults[qIndex * k + i].data << endl;
            }
        }
        file.close();
    }

}
#endif //FILEUTIL_H
