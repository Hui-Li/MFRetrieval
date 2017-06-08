#ifndef INTMATRIXROW_H
#define INTMATRIXROW_H

#include "../util/Base.h"
#include "ExtendMatrixRow.h"

class IntMatrixRow : public ExtendMatrixRow {

public:

    int *iRawData;
    int sumOfCoordinate; //sumOfCoordinate + dimension * 1
    int sumOfCoordinateLeft;
    int sumOfCoordinateRight;

    inline IntMatrixRow(){
        this->iRawData = NULL;
    }

    inline IntMatrixRow(int gRowID, double norm, int colNum):ExtendMatrixRow(gRowID, norm, colNum) {
        this->iRawData = new int[colNum];
    }

    inline IntMatrixRow(int gRowID, int colNum):ExtendMatrixRow(gRowID, colNum) {
        this->iRawData = new int[colNum];
    }


    inline ~IntMatrixRow(){
        if(iRawData) {
            delete[] iRawData;
        }
    }
};


#endif //INTMATRIXROW_H
