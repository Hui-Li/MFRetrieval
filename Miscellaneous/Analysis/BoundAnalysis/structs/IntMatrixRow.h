#ifndef INTMATRIXROW_H
#define INTMATRIXROW_H

#include "../util/Base.h"
#include "ExtendMatrixRow.h"

class IntMatrixRow : public ExtendMatrixRow {

public:

    int *iRawData; // [-128, 127]
#ifdef WITH_SIMD
    int8_t *iRawDataInt8;
#endif

    double intNorm;
    double subIntNorm;
    int sumOfCoordinate; //sumOfCoordinate + dimension * 1

    inline IntMatrixRow(){
        this->iRawData = NULL;

#ifdef WITH_SIMD
        this->iRawDataInt8 = NULL;
#endif

    }

    inline IntMatrixRow(int gRowID, double norm, int colNum):ExtendMatrixRow(gRowID, norm, colNum) {
        this->iRawData = new int[colNum];

#ifdef WITH_SIMD
        this->iRawDataInt8 = new int8_t[colNum];
#endif
    }

    inline ~IntMatrixRow(){
        if(iRawData) {
            delete[] iRawData;
        }

#ifdef WITH_SIMD
        if(iRawDataInt8) {
            delete[] iRawDataInt8;
        }
#endif
    }
};


#endif //INTMATRIXROW_H
