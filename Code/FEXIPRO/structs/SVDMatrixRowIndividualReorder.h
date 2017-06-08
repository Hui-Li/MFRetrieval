#ifndef SVDMATRIXROWINDIVIDUALREORDER_H
#define SVDMATRIXROWINDIVIDUALREORDER_H

#include "../util/Base.h"
#include "ExtendMatrixRow.h"

class SVDMatrixRowIndividualReorder : public ExtendMatrixRow {

public:
    double *vRawData;
//    double vNorm;
    double *subVNorm;

    inline SVDMatrixRowIndividualReorder(){
        this->vRawData = NULL;
        this->subVNorm = NULL;
    }

    inline SVDMatrixRowIndividualReorder(int gRowID, double norm, int colNum):ExtendMatrixRow(gRowID, norm, colNum) {
        this->vRawData = new double[colNum];
        this->subVNorm = new double[colNum];
    }

    inline ~SVDMatrixRowIndividualReorder(){
        if(vRawData) {
            delete[] vRawData;
        }

        if(subVNorm) {
            delete[] subVNorm;
        }
    }
};


#endif //SVDMATRIXROWINDIVIDUALREORDER_H
