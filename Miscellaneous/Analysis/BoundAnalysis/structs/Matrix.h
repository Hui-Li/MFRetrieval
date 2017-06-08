#ifndef MATRIX_H
#define MATRIX_H

#include "../util/Base.h"
#include "VectorElement.h"

class Matrix {

public:
    double *rawData;
    int rowNum;
    int colNum;

    inline Matrix() {
        this->rawData = NULL;
        this->rowNum = 0;
        this->colNum = 0;
    }

    inline ~Matrix() {
        if(rawData){
            delete[] rawData;
        }
    }

    inline double *getRowPtr(const int rowIndex) const {
        return &rawData[rowIndex * colNum];
    }

    inline void init(double *rawData, const int rowNum, const int colNum) {
        this->rowNum = rowNum;
        this->colNum = colNum;
        this->rawData = rawData;
    }

    inline void init(const int rowNum, const int colNum) {
        this->rowNum = rowNum;
        this->colNum = colNum;
        this->rawData = new double[rowNum * colNum];
    }

    inline void setValue(const int rowID, const int colID, const double value){
        this->rawData[rowID * colNum + colID] = value;
    }

    void makePPositive(Matrix &origin, const double maxNorm, vector<VectorElement> &unsortedNorm) {
        this->rowNum = origin.rowNum;
        this->colNum = origin.colNum + 2;
        this->rawData = new double[rowNum * colNum];

        double powMaxNorm = maxNorm * maxNorm;

        // map to d+2
        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {

            double *rowPtr = &rawData[rowIndex * colNum];
            double *oriPtr = origin.getRowPtr(rowIndex);
            double originNorm = unsortedNorm[rowIndex].data;

            rowPtr[1] = sqrt((powMaxNorm - originNorm * originNorm)/powMaxNorm) + 1;

            rowPtr[0] = rowPtr[1] * rowPtr[1];
            for (int colIndex = 2; colIndex < colNum; colIndex++) {
                rowPtr[colIndex] = oriPtr[colIndex-2] / maxNorm + 1;
                rowPtr[0] += rowPtr[colIndex] * rowPtr[colIndex];
            }

        }
    }


    void makePPositive2(Matrix &origin, const double maxNorm, vector<VectorElement> &originalNorms, vector<VectorElement> &newSortedNorms) {
        this->rowNum = origin.rowNum;
        this->colNum = origin.colNum + 2;
        this->rawData = new double[rowNum * colNum];

        double powMaxNorm = maxNorm * maxNorm;

        newSortedNorms.resize(origin.rowNum);

        // map to d+2
        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {

            double *rowPtr = &rawData[rowIndex * colNum];
            double *oriPtr = origin.getRowPtr(rowIndex);

            double originalNorm = originalNorms[rowIndex].data;

            rowPtr[1] = sqrt((powMaxNorm - originalNorm * originalNorm)/powMaxNorm) + 1;

            rowPtr[0] = rowPtr[1] * rowPtr[1];
            for (int colIndex = 2; colIndex < colNum; colIndex++) {
                rowPtr[colIndex] = oriPtr[colIndex-2] / maxNorm + 1;
                rowPtr[0] += rowPtr[colIndex] * rowPtr[colIndex];
            }

            newSortedNorms[rowIndex] = VectorElement(rowIndex, 0.5 * (rowPtr[0] + sqrt(rowPtr[0])));

        }

        sort(newSortedNorms.begin(), newSortedNorms.end(), greater<VectorElement>());
    }

    void readData(string dataFilePath) {
        vector <string> lines;
        string line;

        ifstream fin(dataFilePath.c_str());

        int rowNum = 0;
        int colNum = 0;

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

        this->rowNum = rowNum;
        this->colNum = colNum;
        this->rawData = new double[rowNum * colNum];
        int colIndex = 0;

        for (int rowIndex = 0; rowIndex < lines.size(); rowIndex++) {
            stringstream ss(lines[rowIndex]);
            colIndex = 0;
            while (ss >> this->rawData[rowIndex * colNum + colIndex]) {
                if (ss.peek() == ',') {
                    ss.ignore();
                    colIndex++;
                }
            }
        }
    }
};

#endif //MATRIX_H