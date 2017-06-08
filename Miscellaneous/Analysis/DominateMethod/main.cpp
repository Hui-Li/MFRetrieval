#include <iostream>
#include <set>
#include "Matrix.h"
#include "Calculator.h"

using namespace std;

int main() {
    Matrix p;

    p.readData("p.txt");

    cout << p.rowNum << "," << p.colNum << endl;

    vector<VectorElement> pNorms(p.rowNum);
    double maxNorm = Calculator::calNorms(p, pNorms);

    Matrix positiveP;
    positiveP.makePPositive(p, maxNorm, pNorms);

    vector<int> dominateIDs;

    for (int pID1 = 0; pID1 < positiveP.rowNum; pID1++) {
        const double *pPtr1 = positiveP.getRowPtr(pID1);
        for (int pID2 = 0; pID2 < positiveP.rowNum; pID2++) {
            if (pID2 == pID1) {
                continue;
            }
            const double *pPtr2 = positiveP.getRowPtr(pID2);

            for (int colIndex = 0; colIndex < positiveP.colNum; colIndex++) {
                if (pPtr1[colIndex] < pPtr2[colIndex]) {
                    break;
                }
                if (colIndex == positiveP.colNum - 1) {
                    dominateIDs.push_back(pID2);
                }
            }
        }
    }

    set<int> toDel(dominateIDs.begin(), dominateIDs.end());

    cout << toDel.size() << endl;
    cout << p.rowNum - toDel.size() << endl;

    return 0;
}