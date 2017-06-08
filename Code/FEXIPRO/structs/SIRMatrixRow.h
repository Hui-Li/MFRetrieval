#ifndef SIRMATRIXROW_H
#define SIRMATRIXROW_H

#include "../util/Base.h"
#include "ExtendMatrixRow.h"

class SIRMatrixRow : public ExtendMatrixRow {

public:
    int *iRawData;
    double subNorm;
    int sumOfCoordinate1; //sumOfCoordinate + dimension * 1
    int sumOfCoordinate2;
    double subVNorm;
    double subTransformedSubVNorm;
    double partialSumOfCoordinate;
    double leftPartialSumOfCoordinate;
    double sumOfCoordinate;

    inline SIRMatrixRow(){
        this->iRawData = NULL;
    }

    inline SIRMatrixRow(int gRowID, int colNum, double norm):ExtendMatrixRow(gRowID, norm, colNum) {
        this->iRawData = new int[colNum];
    }

    inline ~SIRMatrixRow(){

        if(iRawData) {
            delete[] iRawData;
        }
    }
};
#endif //SIRMATRIXROW_H