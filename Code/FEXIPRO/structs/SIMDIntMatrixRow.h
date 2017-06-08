#ifndef SIMDINTMATRIXROW_H
#define SIMDINTMATRIXROW_H

#include "../util/Base.h"
#include "ExtendMatrixRow.h"

class SIMDIntMatrixRow : public ExtendMatrixRow {

public:

    int16_t *iRawData;
    int16_t sumOfCoordinate; //sumOfCoordinate + dimension * 1

    inline SIMDIntMatrixRow(){
        this->iRawData = NULL;
    }

    inline SIMDIntMatrixRow(int gRowID, double norm, int colNum):ExtendMatrixRow(gRowID, norm, colNum) {
        this->iRawData = new int16_t[colNum];
    }

    inline ~SIMDIntMatrixRow(){
        if(iRawData) {
            delete[] iRawData;
        }
    }
};


#endif //SIMDINTMATRIXROW_H
