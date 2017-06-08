#ifndef REDSVDMATRIXROW2_H
#define REDSVDMATRIXROW2_H

#include "../util/Base.h"
#include "ExtendMatrixRow.h"

class RedSVDMatrixRow2 : public ExtendMatrixRow {

public:
    double subVNorm;
    double subTransformedSubVNorm;
    double partialSumOfCoordinate;
    double leftPartialSumOfCoordinate;
    double sumOfCoordinate;

    inline RedSVDMatrixRow2(){
    }

    inline RedSVDMatrixRow2(int gRowID, double originalNorm, int colNum):ExtendMatrixRow(gRowID, originalNorm, colNum) {}

};


#endif //REDSVDMATRIXROW2_H
