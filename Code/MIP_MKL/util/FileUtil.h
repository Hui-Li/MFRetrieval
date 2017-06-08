#ifndef FILEUTIL_H
#define FILEUTIL_H

#include "Base.h"
#include "../structs/VectorElement.h"

namespace FileUtil {

    void getSize(string dataFilePath, int &rowNum, int &colNum) {
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

    }

    void readData(string dataFilePath, double *rawData, const int colNum) {

        vector <string> lines;
        string line;

        ifstream fin(dataFilePath.c_str());

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

    void readDataT(string dataFilePath, double *rawData, const int colNum) {

        vector <string> lines;
        string line;

        ifstream fin(dataFilePath.c_str());

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

        int colIndex = 0;

        for (int rowIndex = 0; rowIndex < lines.size(); rowIndex++) {
            stringstream ss(lines[rowIndex]);
            colIndex = 0;

            while (ss >> rawData[colIndex * lines.size() + rowIndex]) {

                if (ss.peek() == ',') {
                    ss.ignore();
                    colIndex++;
                }
            }
        }

    }


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

    void outputResult(vector<vector<VectorElement> > &results, string &resultPath) {
        createFolder(resultPath);
        ofstream file(resultPath.c_str());

        for (int rowIndex = 0; rowIndex < results.size(); rowIndex++) {
            vector<VectorElement> &result = results[rowIndex];
            for (int colIndex = 0; colIndex < result.size(); colIndex++) {
                file << rowIndex << " " << result[colIndex].id << " " << result[colIndex].data << endl;
            }
        }

        file.close();
    }

}
#endif //FILEUTIL_H
