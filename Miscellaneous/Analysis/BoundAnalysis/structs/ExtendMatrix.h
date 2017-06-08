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

    inline void initExtendMatrixSVDIncrPositiveAvg(const Matrix &rawData, const vector<VectorElement> &sortedNorm) {

        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {
            int rowID = sortedNorm[rowIndex].id;
            const double *rawRowPtr = rawData.getRowPtr(rowID);
            this->rows[rowIndex] = *(new T(rowID, sortedNorm[rowID].data, colNum));
            this->rows[rowIndex].rawData[0] = rawRowPtr[0];
            this->rows[rowIndex].rawData[1] = rawRowPtr[1];
            double max = -DBL_MAX;
            for (int colIndex = 2; colIndex < colNum; colIndex++) {
                if(rawRowPtr[colIndex] > max){
                    max = rawRowPtr[colIndex];
                }
                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
            }
            this->rows[rowIndex].l1Norm = max;
        }
    }


    inline void initExtendMatrixWithTwoNorm1(const Matrix &rawData, const vector<VectorElement> &sortedL1Norm, const
    vector<VectorElement> &unsortedL2Norm) {

        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {
            int rowID = sortedL1Norm[rowIndex].id;
            const double *rawRowPtr = rawData.getRowPtr(rowID);
            this->rows[rowIndex] = *(new T(rowID, unsortedL2Norm[rowID].data, sortedL1Norm[rowIndex].data, colNum));

            for (int colIndex = 0; colIndex < colNum; colIndex++) {
                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
            }
        }
    }

    inline void initExtendMatrixWithTwoNorm2(const Matrix &rawData, const vector<VectorElement> &sortedL2Norm, const
    vector<VectorElement> &unsortedL1Norm) {

        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {
            int rowID = sortedL2Norm[rowIndex].id;
            const double *rawRowPtr = rawData.getRowPtr(rowID);
            this->rows[rowIndex] = *(new T(rowID, sortedL2Norm[rowIndex].data, unsortedL1Norm[rowID].data, colNum));

            for (int colIndex = 0; colIndex < colNum; colIndex++) {
                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
            }
        }
    }

    // for svd
    inline void initSVDExtendMatrix(const Matrix &rawData, const Matrix &vRawData, const vector<VectorElement> &sortedNorm,
                        const vector<double> &vNorms, const vector<double> &subVNorms) {
        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {
            const double *rawRowPtr = rawData.getRowPtr(sortedNorm[rowIndex].id);
            const double *vRawRowPtr = vRawData.getRowPtr(sortedNorm[rowIndex].id);
            this->rows[rowIndex] = *(new T(sortedNorm[rowIndex].id, sortedNorm[rowIndex].data, colNum));
            this->rows[rowIndex].vNorm = vNorms[sortedNorm[rowIndex].id];
            this->rows[rowIndex].subVNorm = subVNorms[sortedNorm[rowIndex].id];
            for (int colIndex = 0; colIndex < colNum; colIndex++) {
                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
                this->rows[rowIndex].vRawData[colIndex] = vRawRowPtr[colIndex];
            }

        }
    }

    // for svd + positive
    inline void initSVDExtendMatrix2(const Matrix &rawData, const vector<VectorElement> &sortedNorm, const vector<double> &subNorms) {
        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {
            const double *rawRowPtr = rawData.getRowPtr(sortedNorm[rowIndex].id);
            this->rows[rowIndex] = *(new T(sortedNorm[rowIndex].id, sortedNorm[rowIndex].data, colNum));
            this->rows[rowIndex].l1Norm = subNorms[sortedNorm[rowIndex].id];
            for (int colIndex = 0; colIndex < colNum; colIndex++) {
                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
            }

        }
    }

//    inline void initSVDPositiveExtendMatrix(const Matrix &rawData, const Matrix &vRawData, const vector<VectorElement> &sortedNorm,
//                                    const vector<double> &vNorms, const vector<double> &subVNorms) {
//        this->rowNum = rawData.rowNum;
//        this->colNum = rawData.colNum;
//        this->rows = new T[rowNum];
//
//        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {
//            const double *rawRowPtr = rawData.getRowPtr(sortedNorm[rowIndex].id);
//            const double *vRawRowPtr = vRawData.getRowPtr(sortedNorm[rowIndex].id);
//            this->rows[rowIndex] = *(new T(sortedNorm[rowIndex].id, sortedNorm[rowIndex].data, colNum));
//            this->rows[rowIndex].vNorm = vNorms[sortedNorm[rowIndex].id];
//            this->rows[rowIndex].subVNorm = subVNorms[sortedNorm[rowIndex].id];
//            for (int colIndex = 0; colIndex < colNum; colIndex++) {
//                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
//                this->rows[rowIndex].vRawData[colIndex] = vRawRowPtr[colIndex];
//            }
//
//        }
//    }

    // for int map
    inline void initIntExtendMatrix(const Matrix &rawData, const vector<VectorElement> &sortedNorm, const int mapMaxValue) {
        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        minValue = DBL_MAX;
        maxValue = -DBL_MAX;

        for (int rowID = 0; rowID < rawData.rowNum; rowID++) {

            const double *row = rawData.getRowPtr(rowID);

            for (int colIndex = 0; colIndex < rawData.colNum; colIndex++) {

                if(row[colIndex] < minValue) {
                    minValue = row[colIndex];
                }

                if(row[colIndex] > maxValue) {
                    maxValue = row[colIndex];
                }
            }
        }

        double absMin = fabs(minValue);
        double absMax = fabs(maxValue); // maxValue can be negative
        double denominator = absMin > absMax ? absMin : absMax;
        if(denominator==0){
            ratio = 0;
        } else {
            ratio = mapMaxValue / denominator;
        }

        int sumOfCoordinate; //sumOfCoordinate + dimension * 1
        int temp;
//        double intNorm;
        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {
            const double *rawRowPtr = rawData.getRowPtr(sortedNorm[rowIndex].id);
            sumOfCoordinate = 0;
//            intNorm = 0;
            this->rows[rowIndex] = *(new T(sortedNorm[rowIndex].id, sortedNorm[rowIndex].data, colNum));

            for (int colIndex = 0; colIndex < colNum; colIndex++) {
                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
                temp = floor(rawRowPtr[colIndex] * ratio);
//                intNorm += temp * temp;
                this->rows[rowIndex].iRawData[colIndex] = temp;
#ifdef WITH_SIMD
                this->rows[rowIndex].iRawDataInt8[colIndex] = temp;
#endif
                sumOfCoordinate += fabs(temp);
            }

//            this->rows[rowIndex].intNorm = sqrt(intNorm);
            this->rows[rowIndex].sumOfCoordinate = sumOfCoordinate + colNum;
        }
    }

    // for int + svd
    inline void initIntSVDExtendMatrix(const Matrix &rawData, const Matrix &vRawData, const vector<VectorElement> &sortedNorm,
                                       const int checkDim, const int mapMaxValue) {
        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        minValue = DBL_MAX;
        maxValue = -DBL_MAX;

        int gRowID;
        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {
            gRowID = sortedNorm[rowIndex].id;
            const double *rawRowPtr = rawData.getRowPtr(gRowID);
            const double *vRawRowPtr = vRawData.getRowPtr(gRowID);
            this->rows[rowIndex] = *(new T(gRowID, sortedNorm[rowIndex].data, colNum, colNum));

            double subVNorms = 0;

            for (int colIndex = colNum - 1; colIndex >= checkDim; colIndex--) {

                this->rows[rowIndex].vRawData[colIndex] = vRawRowPtr[colIndex];
                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
                subVNorms += vRawRowPtr[colIndex] * vRawRowPtr[colIndex];

                if(rawRowPtr[colIndex] < minValue) {
                    minValue = rawRowPtr[colIndex];
                }

                if(rawRowPtr[colIndex] > maxValue) {
                    maxValue = rawRowPtr[colIndex];
                }

            }

            this->rows[rowIndex].subVNorm = sqrt(subVNorms);

            for (int colIndex = checkDim - 1; colIndex >= 0; colIndex--) {

                this->rows[rowIndex].vRawData[colIndex] = vRawRowPtr[colIndex];
                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];

                if(rawRowPtr[colIndex] < minValue) {
                    minValue = rawRowPtr[colIndex];
                }

                if(rawRowPtr[colIndex] > maxValue) {
                    maxValue = rawRowPtr[colIndex];
                }

            }

        }

        double absMin = fabs(minValue);
        double absMax = fabs(maxValue); // maxValue can be negative
        double denominator = absMin > absMax ? absMin : absMax;
        if(denominator==0){
            ratio = 0;
        } else {
            ratio = mapMaxValue / denominator;
        }

        int sumOfCoordinate; //sumOfCoordinate + dimension * 1
        int temp;

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {
            const double *rawRowPtr = rawData.getRowPtr(sortedNorm[rowIndex].id);
            sumOfCoordinate = 0;

            for (int colIndex = 0; colIndex < colNum; colIndex++) {
                temp = floor(rawRowPtr[colIndex] * ratio);
                this->rows[rowIndex].iRawData[colIndex] = temp;
                sumOfCoordinate += fabs(temp);
            }

            this->rows[rowIndex].sumOfCoordinate = sumOfCoordinate + colNum;
        }
    }

    // for int map + incr
    inline void initIntExtendMatrix2(const Matrix &rawData, const vector<VectorElement> &sortedNorm, const int mapMaxValue, const int checkDim) {
        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        minValue = DBL_MAX;
        maxValue = -DBL_MAX;

        for (int rowID = 0; rowID < rawData.rowNum; rowID++) {

            const double *row = rawData.getRowPtr(rowID);

            for (int colIndex = 0; colIndex < rawData.colNum; colIndex++) {

                if(row[colIndex] < minValue) {
                    minValue = row[colIndex];
                }

                if(row[colIndex] > maxValue) {
                    maxValue = row[colIndex];
                }
            }
        }

        double absMin = fabs(minValue);
        double absMax = fabs(maxValue); // maxValue can be negative
        double denominator = absMin > absMax ? absMin : absMax;
        if(denominator==0){
            ratio = 0;
        } else {
            ratio = mapMaxValue / denominator;
        }

        int sumOfCoordinate; //sumOfCoordinate + dimension * 1
        int temp;
        double intNorm;
        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {
            const double *rawRowPtr = rawData.getRowPtr(sortedNorm[rowIndex].id);
            sumOfCoordinate = 0;
            intNorm = 0;
            this->rows[rowIndex] = *(new T(sortedNorm[rowIndex].id, sortedNorm[rowIndex].data, colNum));

            for (int colIndex = colNum - 1; colIndex >= 0; colIndex--) {
                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
                temp = floor(rawRowPtr[colIndex] * ratio);
                intNorm += temp * temp;
                this->rows[rowIndex].iRawData[colIndex] = temp;
#ifdef WITH_SIMD
                this->rows[rowIndex].iRawDataInt8[colIndex] = temp;
#endif
                if(colIndex==checkDim) {
                    this->rows[rowIndex].subIntNorm = sqrt(intNorm);
                }

                sumOfCoordinate += fabs(temp);
            }

//            this->rows[rowIndex].intNorm = sqrt(intNorm);
            this->rows[rowIndex].sumOfCoordinate = sumOfCoordinate + colNum;
        }
    }

    // for svd int map
    inline void initSVDIntExtendMatrix(const Matrix &rawData, const Matrix &vRawData,
                                       const vector<VectorElement> &sortedNorm,
                                       const int checkDim) {
        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        double vNorm = 0;
        double subVNorm = 0;
        int currentGID;
        minValue = DBL_MAX;
        maxValue = -DBL_MAX;

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {

            currentGID = sortedNorm[rowIndex].id;
            const double *rawRowPtr = rawData.getRowPtr(currentGID);
            const double *vRawRowPtr = vRawData.getRowPtr(currentGID);

            this->rows[rowIndex] = *(new T(sortedNorm[rowIndex].id, sortedNorm[rowIndex].data, colNum));

            for (int colIndex = colNum - 1; colIndex >= 0; colIndex--) {

                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
                this->rows[rowIndex].vRawData[colIndex] = vRawRowPtr[colIndex];

                vNorm += vRawRowPtr[colIndex] * vRawRowPtr[colIndex];

                if (vRawRowPtr[colIndex] < minValue) {
                    minValue = vRawRowPtr[colIndex];
                }

                if (vRawRowPtr[colIndex] > maxValue) {
                    maxValue = vRawRowPtr[colIndex];
                }

                if(colIndex==checkDim) {
                    subVNorm = sqrt(vNorm);
                }
            }

            this->rows[rowIndex].vNorm = sqrt(vNorm);
            this->rows[rowIndex].subVNorm = subVNorm;
        }

        double absMin = fabs(minValue);
        double absMax = fabs(maxValue); // maxValue can be negative
        double denominator = absMin > absMax ? absMin : absMax;
        if(denominator==0){
            ratio = 0;
        } else {
            // ToDo: hard code!!!
            ratio = 127 / denominator;
        }

        int sumOfCoordinate; //sumOfCoordinate + dimension * 1
        int temp;

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {

            currentGID = sortedNorm[rowIndex].id;
            const double *vRawRowPtr = vRawData.getRowPtr(currentGID);

            sumOfCoordinate = 0;

            for (int colIndex = colNum - 1; colIndex >= checkDim; colIndex--) {

                temp = floor(vRawRowPtr[colIndex] * ratio);

                this->rows[rowIndex].iRawData[colIndex] = temp;

                sumOfCoordinate += fabs(temp);

                if(colIndex==checkDim) {
                    break;
                }
            }

            sumOfCoordinate += colNum - checkDim;

            this->rows[rowIndex].sumOfCoordinate = sumOfCoordinate;

        }

    }

    // for svd int map in first part
    inline void initSVDIntExtendMatrixNew(const Matrix &rawData, const Matrix &vRawData,
                                       const vector<VectorElement> &sortedNorm,
                                       const int checkDim) {
        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        double vNorm = 0;
        double subVNorm = 0;
        int currentGID;

        minValue = DBL_MAX;
        maxValue = -DBL_MAX;

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {

            currentGID = sortedNorm[rowIndex].id;
            const double *rawRowPtr = rawData.getRowPtr(currentGID);
            const double *vRawRowPtr = vRawData.getRowPtr(currentGID);

            this->rows[rowIndex] = *(new T(sortedNorm[rowIndex].id, sortedNorm[rowIndex].data, colNum, checkDim));

            for (int colIndex = colNum - 1; colIndex >= checkDim; colIndex--) {

                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
                this->rows[rowIndex].vRawData[colIndex] = vRawRowPtr[colIndex];

                vNorm += vRawRowPtr[colIndex] * vRawRowPtr[colIndex];

            }

            subVNorm = sqrt(vNorm);

//            sumOfCoordinate = 0;

            for (int colIndex = checkDim-1; colIndex >= 0; colIndex--) {
                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
                this->rows[rowIndex].vRawData[colIndex] = vRawRowPtr[colIndex];

                vNorm += vRawRowPtr[colIndex] * vRawRowPtr[colIndex];

                if (vRawRowPtr[colIndex] < minValue) {
                    minValue = vRawRowPtr[colIndex];
                }

                if (vRawRowPtr[colIndex] > maxValue) {
                    maxValue = vRawRowPtr[colIndex];
                }

            }

            this->rows[rowIndex].vNorm = sqrt(vNorm);
            this->rows[rowIndex].subVNorm = subVNorm;

        }

        double absMin = fabs(minValue);
        double absMax = fabs(maxValue); // maxValue can be negative
        double denominator = absMin > absMax ? absMin : absMax;
        if(denominator==0){
            ratio = 0;
        } else {
            // ToDo: hard code!!!
            ratio = 127 / denominator;
        }

        int sumOfCoordinate; //sumOfCoordinate + dimension * 1
        int temp;

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {

            currentGID = sortedNorm[rowIndex].id;
            const double *vRawRowPtr = vRawData.getRowPtr(currentGID);

            sumOfCoordinate = 0;

            for (int colIndex = 0; colIndex < checkDim; colIndex++) {

                temp = floor(vRawRowPtr[colIndex] * ratio);

                this->rows[rowIndex].iRawData[colIndex] = temp;

                sumOfCoordinate += fabs(temp);

            }

            sumOfCoordinate += checkDim;

            this->rows[rowIndex].sumOfCoordinate = sumOfCoordinate;

        }

    }

    // for svd int all
    inline void initSVDIntExtendMatrix2(const Matrix &rawData, const Matrix &vRawData,
                                       const vector<VectorElement> &sortedNorm,
                                       const int checkDim) {
        this->rowNum = rawData.rowNum;
        this->colNum = rawData.colNum;
        this->rows = new T[rowNum];

        double vNorm = 0;
        double subVNorm = 0;
        int currentGID;
        minValue = DBL_MAX;
        maxValue = -DBL_MAX;

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {

            currentGID = sortedNorm[rowIndex].id;
            const double *rawRowPtr = rawData.getRowPtr(currentGID);
            const double *vRawRowPtr = vRawData.getRowPtr(currentGID);

            this->rows[rowIndex] = *(new T(sortedNorm[rowIndex].id, sortedNorm[rowIndex].data, colNum));

            for (int colIndex = colNum - 1; colIndex >= 0; colIndex--) {

                this->rows[rowIndex].rawData[colIndex] = rawRowPtr[colIndex];
                this->rows[rowIndex].vRawData[colIndex] = vRawRowPtr[colIndex];

                vNorm += vRawRowPtr[colIndex] * vRawRowPtr[colIndex];

                if (vRawRowPtr[colIndex] < minValue) {
                    minValue = vRawRowPtr[colIndex];
                }

                if (vRawRowPtr[colIndex] > maxValue) {
                    maxValue = vRawRowPtr[colIndex];
                }

                if(colIndex==checkDim) {
                    subVNorm = sqrt(vNorm);
                }
            }

            this->rows[rowIndex].vNorm = sqrt(vNorm);
            this->rows[rowIndex].subVNorm = subVNorm;
        }

        double absMin = fabs(minValue);
        double absMax = fabs(maxValue); // maxValue can be negative
        double denominator = absMin > absMax ? absMin : absMax;
        if(denominator==0){
            ratio = 0;
        } else {
            ratio = 127 / denominator;
        }

        int sumOfCoordinate; //sumOfCoordinate + dimension * 1
        int temp;

        double iNorm = 0;
        double subINorm = 0;

        for (int rowIndex = 0; rowIndex < rowNum; rowIndex++) {

            currentGID = sortedNorm[rowIndex].id;
            const double *vRawRowPtr = vRawData.getRowPtr(currentGID);

            sumOfCoordinate = 0;
            iNorm = 0;
            subINorm = 0;

            for (int colIndex = colNum - 1; colIndex >= 0; colIndex--) {

                temp = floor(vRawRowPtr[colIndex] * ratio);

                this->rows[rowIndex].iRawData[colIndex] = temp;

                iNorm += temp * temp;
                sumOfCoordinate += fabs(temp);

                if(colIndex==checkDim) {
                    subINorm = sqrt(iNorm);
                }
            }

            sumOfCoordinate += colNum;

            this->rows[rowIndex].sumOfCoordinate = sumOfCoordinate;
            this->rows[rowIndex].iNorm = sqrt(iNorm);
            this->rows[rowIndex].subINorm = subINorm;
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
