#ifndef SVDMATRIXROW_H
#define SVDMATRIXROW_H

#include "../util/Base.h"
#include "ExtendMatrixRow.h"

class SVDMatrixRow : public ExtendMatrixRow {

public:
    double *vRawData;
    double vNorm;
    double subVNorm;

    inline SVDMatrixRow(){
        this->vRawData = NULL;
    }

    inline SVDMatrixRow(int gRowID, double norm, int colNum):ExtendMatrixRow(gRowID, norm, colNum) {
        this->vRawData = new double[colNum];
    }

    inline SVDMatrixRow(int gRowID, double originalNorm, int originalColNum, int newColNum):ExtendMatrixRow(gRowID, originalNorm, originalColNum) {
        this->vRawData = new double[newColNum];
    }

    inline ~SVDMatrixRow(){
        if(vRawData) {
            delete[] vRawData;
        }
    }
};


#endif //SVDMATRIXROW_H
