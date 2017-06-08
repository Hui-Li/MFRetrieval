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

    inline double l2Distance(const double *qRow, const double *pRow, const int dim) {

        double value = 0;

        for (int colIndex = 0; colIndex < dim; colIndex++) {
            value += (qRow[colIndex] - pRow[colIndex]) * (qRow[colIndex] - pRow[colIndex]);
        }

        return sqrt(value);
    }

    inline void calSingleNorm(const double *ptr, const int dim, double &norm) {

        norm = 0;
        for (int i = 0; i < dim; i++) {
            norm += ptr[i] * ptr[i];
        }
        norm = sqrt(norm);
    }

    inline void calNorms(const Matrix &m, vector<VectorElement> &norms, double &maxNorm) {
        norms.resize(m.rowNum);
        maxNorm = -1;
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

    }

    inline void calNorms(const Matrix &m, vector<VectorElement> &norms, double &maxNorm, double &minValue) {
        norms.resize(m.rowNum);
        minValue = DBL_MAX;
        maxNorm = -1;
        for (int rowID = 0; rowID < m.rowNum; rowID++) {
            double norm = 0;
            const double *row = m.getRowPtr(rowID);
            for (int colIndex = 0; colIndex < m.colNum; colIndex++) {
                norm += row[colIndex] * row[colIndex];
                if(row[colIndex] < minValue) {
                    minValue = row[colIndex];
                }
            }
            norm = sqrt(norm);
            norms[rowID] = VectorElement(rowID, norm);
            if (norm > maxNorm) {
                maxNorm = norm;
            }
        }


        if(minValue < 0){
            minValue = abs(minValue);
        } else {
            minValue = 0;
        }

    }


    inline void calNorms(const Matrix &m, vector<VectorElement> &norms) {
        norms.resize(m.rowNum);

        for (int rowID = 0; rowID < m.rowNum; rowID++) {
            double norm = 0;
            const double *row = m.getRowPtr(rowID);
            for (int colIndex = 0; colIndex < m.colNum; colIndex++) {
                norm += row[colIndex] * row[colIndex];
            }
            norm = sqrt(norm);
            norms[rowID] = VectorElement(rowID, norm);
        }

    }

}
#endif //CALCULATOR_H
