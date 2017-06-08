#ifndef TRANSFORMUTIL_H
#define TRANSFORMUTIL_H

#include "Base.h"
#include "../structs/VectorElement.h"
#include "../structs/Matrix.h"

namespace TransformUtil {

    void makePPositive(const Matrix &origin, Matrix &newMatrix, const double maxNorm, const int checkDim,
                       const vector <VectorElement> &unsortedOriginalNorm, vector <VectorElement> &sortedNewNorm,
                       vector<double> &unsortedSubNorm, const vector<double> &addend) {

        newMatrix.rowNum = origin.rowNum;
        newMatrix.colNum = origin.colNum + 2;
        newMatrix.rawData = new double[newMatrix.rowNum * newMatrix.colNum];

        double powMaxNorm = maxNorm * maxNorm;
        sortedNewNorm.resize(origin.rowNum);
        unsortedSubNorm.resize(origin.rowNum);

        double temp;

        for (int rowIndex = 0; rowIndex < newMatrix.rowNum; rowIndex++) {
            double *rowPtr = newMatrix.getRowPtr(rowIndex);
            double *oriPtr = origin.getRowPtr(rowIndex);
            double originNorm = unsortedOriginalNorm[rowIndex].data;

            rowPtr[1] = sqrt((powMaxNorm - originNorm * originNorm)/powMaxNorm) + 1;
            rowPtr[0] = rowPtr[1] * rowPtr[1];

            double norm = 0;

            for (int colIndex = newMatrix.colNum - 1; colIndex >= 2; colIndex--) {
                rowPtr[colIndex] = oriPtr[colIndex-2] / maxNorm + addend[colIndex-2];
                temp = rowPtr[colIndex] * rowPtr[colIndex];
                rowPtr[0] += temp;
                norm += temp;

                if(colIndex==checkDim){
                    unsortedSubNorm[rowIndex] = sqrt(norm);
                }

            }

            sortedNewNorm[rowIndex] = VectorElement(rowIndex, 0.5 * (norm + rowPtr[1] * rowPtr[1] + rowPtr[0] * rowPtr[0])/(norm + rowPtr[1] * rowPtr[1]));
        }

        sort(sortedNewNorm.begin(), sortedNewNorm.end(), greater<VectorElement>());
    }

//    void makePPositiveSVDIncr(const Matrix &vRawData, Matrix &newMatrix, const double maxVNorm, const int checkDim,
//                       const vector <VectorElement> &sortedOriginalNorm, const vector <VectorElement> &unsortedVNorms,
//                       vector<double> &sortedSubNorm, const vector<double> &addend) {
//
//        newMatrix.rowNum = vRawData.rowNum;
//        newMatrix.colNum = vRawData.colNum + 2;
//        newMatrix.rawData = new double[newMatrix.rowNum * newMatrix.colNum];
//
//        double powMaxNorm = maxVNorm * maxVNorm;
//        sortedSubNorm.resize(vRawData.rowNum);
//
//        double temp;
//
//        for (int rowIndex = 0; rowIndex < newMatrix.rowNum; rowIndex++) {
//            int gRowID = sortedOriginalNorm[rowIndex].id;
//            double *rowPtr = newMatrix.getRowPtr(rowIndex);
//            double *oriPtr = vRawData.getRowPtr(gRowID);
//            double originVNorm = unsortedVNorms[gRowID].data;
//
//            rowPtr[1] = sqrt((powMaxNorm - originVNorm * originVNorm)/powMaxNorm) + 1;
//            rowPtr[0] = rowPtr[1] * rowPtr[1];
//
//            double norm = 0;
//
//            for (int colIndex = newMatrix.colNum - 1; colIndex >= 2; colIndex--) {
//                rowPtr[colIndex] = oriPtr[colIndex-2] / maxVNorm + addend[colIndex-2];
//                temp = rowPtr[colIndex] * rowPtr[colIndex];
//                rowPtr[0] += temp;
//                norm += temp;
//
//                if(colIndex==checkDim){
//                    sortedSubNorm[rowIndex] = sqrt(norm);
//                }
//
//            }
//
//        }
//
//    }

    void makePPositiveSVDIncr(const Matrix &vRawData, Matrix &newMatrix, const double maxVNorm, const int checkDim,
                              const vector <VectorElement> &sortedOriginalNorm, const vector <VectorElement> &unsortedVNorms,
                              vector<double> &sortedSubNorm, vector<double> &addend) {

        newMatrix.rowNum = vRawData.rowNum;
        newMatrix.colNum = vRawData.colNum + 2;
        newMatrix.rawData = new double[newMatrix.rowNum * newMatrix.colNum];

        double powMaxNorm = maxVNorm * maxVNorm;
        sortedSubNorm.resize(vRawData.rowNum);

        for (int rowIndex = 0; rowIndex < newMatrix.rowNum; rowIndex++) {
            int gRowID = sortedOriginalNorm[rowIndex].id;
            double *rowPtr = newMatrix.getRowPtr(rowIndex);
            double *oriPtr = vRawData.getRowPtr(gRowID);
            double originVNorm = unsortedVNorms[gRowID].data;

            rowPtr[1] = sqrt(powMaxNorm - originVNorm * originVNorm);

            double norm = 0;

            for (int colIndex = newMatrix.colNum - 1; colIndex >= 2; colIndex--) {
                rowPtr[colIndex] = oriPtr[colIndex-2] + addend[colIndex-2];
                norm += rowPtr[colIndex] * rowPtr[colIndex];

                if(colIndex==checkDim){
                    sortedSubNorm[rowIndex] = sqrt(norm);
                }

            }

            rowPtr[0] = rowPtr[1] * rowPtr[1] + norm;
        }

    }
}

#endif //TRANSFORMUTIL_H
