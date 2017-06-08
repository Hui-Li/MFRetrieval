#include "util/Base.h"
#include "IntScale.h"

using namespace std;

int main() {

//    string dataSets[4] = {"MovieLens", "Yelp", "Netflix", "KDD"};
    string dataSets[3] = {"MovieLens", "Yelp", "Netflix"};
//    string dataSets[1] = {"Netflix"};

    for(auto data:dataSets){

        cout << "---------------------" << endl;
        cout << data << endl;

        // for max dimension
        IntScale::analysis(data);
    }

    return 0;
}