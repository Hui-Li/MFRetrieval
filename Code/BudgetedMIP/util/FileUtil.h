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

    void outputResult(const int k, const int d, const int qNum, const double *q, const double *p, vector<vector<size_t> > &candidates, string &resultPath) {

        ofstream file(resultPath.c_str());

        for (int qIndex = 0; qIndex < qNum; qIndex++) {
            vector<size_t > &candidate = candidates[qIndex];
            for (int i = 0; i < k; i++) {
                const double *qPtr = q + qIndex * d;
                const double *pPtr = p + candidate[i] * d;
                double ip = std::inner_product(qPtr, qPtr + d, pPtr, 0.0);
                file << qIndex << "," << candidate[i] << "," << ip << endl;
            }
        }
        file.close();
    }

    void readGroundtruth(vector<unordered_set<int> > &groundTruth, string groundTruthPath) {

        string line;

        ifstream fin(groundTruthPath.c_str());

        vector<string> pars;
        while (getline(fin, line)) {
            if (line.length() == 0) {
                continue;
            }

            boost::split(pars, line, boost::is_any_of(","));
            int qID = stoi(pars[0]);
            int pID = stoi(pars[1]);
            groundTruth[qID].insert(pID);
        }

        fin.close();
    }

}
#endif //FILEUTIL_H
