#ifndef CALCULATOR_H
#define CALCULATOR_H

#include "Base.h"
#include "../structs/Matrix.h"
#include "../structs/VectorElement.h"

namespace Calculator{


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

    inline double innerProduct(const double *qRow, const double *pRow, const int dim) {

        double value = 0;

        for (int colIndex = 0; colIndex < dim; colIndex++) {
            value += qRow[colIndex] * pRow[colIndex];
        }

        return value;
    }

}
#endif //CALCULATOR_H
