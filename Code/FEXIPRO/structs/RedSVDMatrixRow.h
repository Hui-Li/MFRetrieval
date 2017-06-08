#ifndef REDSVDMATRIXROW_H
#define REDSVDMATRIXROW_H

#include "../util/Base.h"
#include "ExtendMatrixRow.h"

class RedSVDMatrixRow : public ExtendMatrixRow {

public:
    double subVNorm;
    double partialSumOfCoordinate;
    double sumOfCoordinate;

    inline RedSVDMatrixRow(){
    }

    inline RedSVDMatrixRow(int gRowID, double originalNorm, int colNum):ExtendMatrixRow(gRowID, originalNorm, colNum) {}

};


#endif //REDSVDMATRIXROW_H
