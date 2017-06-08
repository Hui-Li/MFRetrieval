#ifndef BOUNDANALYSIS_INTSCALE_H
#define BOUNDANALYSIS_INTSCALE_H

#include "alg/int/MapIntBasicPrune.h"
#include "util/Base.h"
#include <random>

namespace IntScale {

    void analysis(string data){

        Matrix *p = new Matrix();
        p->readData("../../../data/" + data + "/p.txt");
        Matrix *q = new Matrix();
        q->readData("../../../data/"+ data + "/q.txt");

        cout << "P:" << p->rowNum << "," << p->colNum << endl;
        cout << "Q:" << q->rowNum << "," << q->colNum << endl;

        int k = 1;
        int maxScaleValues[6] = {1, 10, 100, 127, 1000, 10000};
//        int qSampleSize = 1000;
//
//        random_device rd;     // only used once to initialise (seed) engine
//        mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
//        uniform_int_distribution<int> uni(0, q->rowNum-1); // guaranteed unbiased
//
//        vector<int> sampleQIDs(qSampleSize);
//        for (int i = 0; i < qSampleSize; i++) {
//            sampleQIDs[i] = uni(rng);
//        }

        for(auto maxScaleValue:maxScaleValues) {
            cout << "----------------------------" << endl;
            cout << "maxScaleValue: " << maxScaleValue << endl;
            MapIntNormPune mapIntNormPune(k, maxScaleValue, q, p);
            mapIntNormPune.topK();
        }

        delete p;
        delete q;
    }
}
#endif //BOUNDANALYSIS_INTSCALE_H
