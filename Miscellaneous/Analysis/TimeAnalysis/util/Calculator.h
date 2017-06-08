#ifndef CALCULATOR_H
#define CALCULATOR_H

#include "Base.h"
#include "../structs/Matrix.h"
#include "../structs/VectorElement.h"

namespace Calculator{

    inline double innerProduct(const double *qRow, const double *pRow, const int dim) {

        double value = 0;

        for (int colIndex = 0; colIndex < dim; colIndex++) {
            value += qRow[colIndex] * pRow[colIndex];
        }

        return value;
    }

    inline double innerProductForTest(volatile const double *qRow, volatile const double *pRow, const int dim) {

        volatile double value = 0;

        for (int colIndex = 0; colIndex < dim; colIndex++) {
//            value += qRow[colIndex] * pRow[colIndex];
            value = qRow[colIndex];
            value = pRow[colIndex];
        }

        return value;
    }


    inline void calSingleNorm(const double *ptr, const int dim, double &norm) {

        norm = 0;
        for (int i = 0; i < dim; i++) {
            norm += ptr[i] * ptr[i];
        }
        norm = sqrt(norm);
    }

    inline void calSingleNorm(const double *ptr, const int dim, const int checkDim, double &norm, double &subNorm) {

        norm = 0;
        subNorm = 0;
        for (int i = 0; i < dim; i++) {
            norm += ptr[i] * ptr[i];

            if (i == checkDim) {
                subNorm = sqrt(norm);
            }

        }
        norm = sqrt(norm);
    }


    inline void calSingleNormForTest(volatile const double *ptr, volatile const int dim, volatile double &norm) {

        norm = 0;
        for (int i = 0; i < dim; i++) {
            norm = ptr[i];
            norm = ptr[i];
//            norm += ptr[i] * ptr[i];
        }
//        norm = sqrt(norm);
        norm = ptr[0];

    }

    inline double calNorms(const Matrix &m, vector<VectorElement> &norms) {
        norms.resize(m.rowNum);
        double maxNorm = -1;
        for (int rowID = 0; rowID < m.rowNum; rowID++) {
            double norm = 0;
            const double *row = m.getRowPtr(rowID);
            for (int colIndex = 0; colIndex < m.colNum; colIndex++) {
                norm += row[colIndex] * row[colIndex];
            }
            norm = sqrt(norm);
            norms[rowID] = VectorElement(rowID, norm);
            if (norm > maxNorm) {
                maxNorm = norm;
            }
        }
        return maxNorm;
    }

    inline double calNormAndSubNorm(const Matrix &m, const int checkDim, vector<VectorElement> &norms, vector<double> &subNorms) {
        norms.resize(m.rowNum);
        subNorms.resize(m.rowNum);

        double maxNorm = -1;
        for (int rowID = 0; rowID < m.rowNum; rowID++) {
            double norm = 0;
            const double *row = m.getRowPtr(rowID);
            for (int colIndex = 0; colIndex < m.colNum; colIndex++) {
                norm += row[colIndex] * row[colIndex];

                if (colIndex == checkDim) {
                    subNorms[rowID] = sqrt(norm);
                }

            }
            norm = sqrt(norm);
            norms[rowID] = VectorElement(rowID, norm);
            if (norm > maxNorm) {
                maxNorm = norm;
            }
        }

        return maxNorm;
    }

    inline void calNormsForTest(const Matrix &m, vector<VectorElement> &norms) {
        norms.resize(m.rowNum);

        for (int rowID = 0; rowID < m.rowNum; rowID++) {
            double norm = 0;
            volatile double *row = m.getRowPtr(rowID);
            for (int colIndex = 0; colIndex < m.colNum; colIndex++) {
                row[colIndex];
                row[colIndex];
//                norm += row[colIndex] * row[colIndex];
            }
//            norm = sqrt(norm);
            norm;
            norm;
            norms[rowID] = VectorElement(rowID, norm);

        }
    }


}
#endif //CALCULATOR_H
