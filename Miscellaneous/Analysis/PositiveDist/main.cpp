#include <iostream>
#include "Matrix.h"

using namespace std;

inline double calNorms(const Matrix &m, vector<VectorElement> &norms) {
    norms.resize(m.rowNum);
    double maxNorm = -1;
    for (int rowID = 0; rowID < m.rowNum; rowID++) {
        double norm = 0;
        const double *row = m.getRowPtr(rowID);
        for (int colIndex = 0; colIndex < m.colNum; colIndex++) {
            norm += row[colIndex] * row[colIndex];
        }
        norm = sqrt(norm);
        norms[rowID] = VectorElement(rowID, norm);
        if (norm > maxNorm) {
            maxNorm = norm;
        }
    }
    return maxNorm;
}

int main() {
    Matrix p;
    p.readData("p.txt");
    Matrix q;
    q.readData("q.txt");

    cout << p.rowNum << "," << p.colNum << endl;
    cout << q.rowNum << "," << q.colNum << endl;

    vector<VectorElement> pNorms;
    double maxPNorm = calNorms(p, pNorms);

    Matrix positiveP;
    positiveP.makePPositive(p, maxPNorm, pNorms);

    vector<VectorElement> qNorms;
    calNorms(q, qNorms);
    Matrix positiveQ;
    positiveQ.makeQPositive(q, qNorms);

    string pPath = "positveP.txt";
    string qPath = "positveQ.txt";

    ofstream file(pPath.c_str());
    for (int i = 0; i < positiveP.rowNum; i++) {
        const double *rowPtr = positiveP.getRowPtr(i);
        for (int j = 1; j < positiveP.colNum; j++) {

            file << rowPtr[j] << endl;

        }
    }

    file.close();

    ofstream file2(qPath.c_str());
    for (int i = 0; i < positiveQ.rowNum; i++) {
        const double *rowPtr = positiveQ.getRowPtr(i);
        for (int j = 1; j < positiveQ.colNum; j++) {

            file2 << rowPtr[j] << endl;

        }
    }

    file2.close();

    return 0;
}