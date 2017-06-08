#ifndef EXTENDMATRIX_H
#define EXTENDMATRIX_H

#include <bitset>
#include "../util/Base.h"
#include "VectorElement.h"


// sorted matrix by norm
template <class T>
class ExtendMatrix {

private:
    T *rows;

public:
    int rowNum;
    int colNum;
    // for int map
    double minValue;
    double maxValue;
    double ratio;

    inline ExtendMatrix(){
        this->rows = NULL;
    }

    inline void initExtendMatrix(const Matrix &rawData, const vector<VectorElement> &sortedNorm) {
        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {
            const double *rawRowPtr = rawData.getRowPtr(sortedNorm[rowIndex].id);
            this->rows[rowIndex] = *(new T(sortedNorm[rowIndex].id, sortedNorm[rowIndex].data, colNum));
            for (int colIndex = 0; colIndex < colNum; colIndex++) {
                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
            }
        }
    }

    inline void initExtendMatrix(const Matrix &rawData, const vector<VectorElement> &sortedNorm, const vector<double> &subNroms) {
        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {
            const double *rawRowPtr = rawData.getRowPtr(sortedNorm[rowIndex].id);
            this->rows[rowIndex] = *(new T(sortedNorm[rowIndex].id, sortedNorm[rowIndex].data, subNroms[sortedNorm[rowIndex].id], colNum));

            for (int colIndex = 0; colIndex < colNum; colIndex++) {
                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
            }
        }
    }


    inline ~ExtendMatrix(){
        if(rows) {
            delete[] rows;
        }
    }

    inline T *getRowPtr(const int rowIndex) const {
        return &rows[rowIndex];
    }


};

#endif //EXTENDMATRIX_H
