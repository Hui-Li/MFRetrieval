#include <iostream>
#include "SVDDisAnalysis.h"
#include "DisAnalysis.h"
#include "SVDDisAnalysis2.h"

using namespace std;

int main(){

    string dataSets[4] = {"MovieLens", "Yelp", "Netflix", "KDD"};

//    string dataSets[1] = {"MovieLens"};

    for(auto data:dataSets){

        cout << "---------------------" << endl;
        cout << data << endl;

        // for original values
        DisAnalysis::avgVector(data);

        // for reorder values
//        DisAnalysis::avgReorderedVector(data);

        // for max dimension
//        DisAnalysis::maxFrequencyVector(data);

        // for svd
//        SVDDisAnalysis::analysis(data);
//        SVDDisAnalysis::analysisForLength(data);

        // for svd2
//        SVDDisAnalysis2::analysis(data);
    }

    return 0;
}