#ifndef VALUEDISTRIBUTION_SVDDISANALYSIS2_H
#define VALUEDISTRIBUTION_SVDDISANALYSIS2_H

#include <armadillo>
#include "Matrix.h"
using namespace std;
using namespace arma;

namespace SVDDisAnalysis2 {

    void SVD(const mat &P_t, const int m, const int n, Matrix &Q, Matrix &newQ, Matrix &newP, Matrix &uForQ) {

        mat U_t;
        vec s;
        mat V;

        svd_econ(U_t, s, V, P_t);

        U_t = U_t.t();

        double *uDataForQ = new double[m * m];
        double *uDataForP = new double[m * m];
        double *vData = new double[m * n];

        for (int rowIndex = 0; rowIndex < m; rowIndex++) {
            for (int colIndex = 0; colIndex < m; colIndex++) {
                uDataForQ[rowIndex * m + colIndex] = U_t(rowIndex, colIndex);
            }
        }

        for (int rowIndex = 0; rowIndex < m; rowIndex++) {
            for (int colIndex = 0; colIndex < m; colIndex++) {
                uDataForP[rowIndex * m + colIndex] = s[rowIndex] * U_t(rowIndex, colIndex);
            }
        }

        uForQ.init(uDataForQ, m, m);

        for (int QIndex = 0; QIndex < newQ.rowNum; QIndex++) {

            double *newQPtr = newQ.getRowPtr(QIndex);
            double *qPtr = Q.getRowPtr(QIndex);
            for (int rowIndex = uForQ.rowNum - 1; rowIndex >= 0; rowIndex--) {
                newQPtr[rowIndex] = 0;
                const double *uPtr = uForQ.getRowPtr(rowIndex);

                for (int colIndex = 0; colIndex < uForQ.colNum; colIndex++) {
                    newQPtr[rowIndex] += uPtr[colIndex] * qPtr[colIndex];
                }
            }
        }

        for(int rowIndex = 0; rowIndex < n; rowIndex++){
            for (int colIndex = 0; colIndex < m; colIndex++) {
                vData[rowIndex * m + colIndex] = V(rowIndex, colIndex);
            }
        }


        double *pData = new double[n * m];

        for(int rowIndex = 0; rowIndex < n; rowIndex++){
            for (int colIndex = 0; colIndex < m; colIndex++) {
                pData[rowIndex * m + colIndex] = 0;
            }
        }

        for (int i = 0; i < n; i++) {
            for (int k = 0; k < m; k++) {
                for (int j = 0; j < m; j++) {
                    pData[i * m + k] += uDataForP[i* m + k] * vData[k * m + j];
                }
            }
        }
        delete[] uDataForP;
        newP.init(pData, n, m);

    }

    void transferQ(const Matrix *u, const double *qPtr, double *newQ) {
        for (int rowIndex = u->rowNum - 1; rowIndex >= 0; rowIndex--) {

            newQ[rowIndex] = 0;
            const double *uPtr = u->getRowPtr(rowIndex);

            for (int colIndex = 0; colIndex < u->colNum; colIndex++) {
                newQ[rowIndex] += uPtr[colIndex] * qPtr[colIndex];
            }

        }

    }

    void analysis(string data){
        mat P_t;
        P_t.load("../../../data/" + data + "/p.txt", csv_ascii);
        P_t = P_t.t();
        cout << P_t.n_rows << "," << P_t.n_cols << endl;

        Matrix *P = new Matrix();
        P->readData("../../../data/" + data + "/p.txt");
        Matrix *Q = new Matrix();
        Q->readData("../../../data/"+ data + "/q.txt");

        Matrix *newQ = new Matrix(Q->rowNum, Q->colNum);
        Matrix *newP = new Matrix();
        Matrix *uForQ = new Matrix();
        SVD(P_t, P->colNum, P->rowNum, *Q, *newQ, *newP, *uForQ);

//        vector<double> avg(newP->colNum, 0);
//        for (int rowIndex = 0; rowIndex < newP->rowNum; rowIndex++) {
//            const double *ptr = newP->getRowPtr(rowIndex);
//            for (int colIndex = 0; colIndex < newP->colNum; colIndex++) {
//                avg[colIndex] += fabs(ptr[colIndex]);
//            }
//        }
//
//        for (int colIndex = 0; colIndex < newP->colNum; colIndex++) {
//            avg[colIndex] /= newP->rowNum;
//        }
//
//        cout << "P:" << newP->rowNum << "," << newP->colNum << endl;
//        string pFile = data + "_pSVDAvg.txt";
//        ofstream outputFile(pFile.c_str());
//
//        for(int colIndex = 0; colIndex < newP->colNum; colIndex++) {
//            outputFile << colIndex + 1 << "," << avg[colIndex] << endl;
//        }
//
//        outputFile.close();
//
//        vector<double> avg2(newQ->colNum, 0);
//        for (int rowIndex = 0; rowIndex < newQ->rowNum; rowIndex++) {
//            const double *qtr = newQ->getRowPtr(rowIndex);
//            for (int colIndex = 0; colIndex < newQ->colNum; colIndex++) {
//                avg2[colIndex] += fabs(qtr[colIndex]);
//            }
//        }
//
//        for (int colIndex = 0; colIndex < newQ->colNum; colIndex++) {
//            avg2[colIndex] /= newQ->rowNum;
//        }
//
//        cout << "Q:" << newQ->rowNum << "," << newQ->colNum << endl;
//        string qFile = data + "_qSVDAvg.txt";
//        ofstream outputFile2(qFile.c_str());
//
//        for(int colIndex = 0; colIndex < newQ->colNum; colIndex++) {
//            outputFile2 << colIndex + 1 << "," << avg2[colIndex] << endl;
//        }
//
//        outputFile2.close();

        const double *qPtr = Q->getRowPtr(131);
        const double *pPtr = P->getRowPtr(444);

        const double *newQPtr = newQ->getRowPtr(131);
        const double *newPPtr = newP->getRowPtr(444);

        double *testQ = new double[Q->colNum];
        transferQ(uForQ, qPtr, testQ);

        double v1 = 0;
        double v2 = 0;
        double v3 = 0;
        for (int i = 0; i < newQ->colNum; i++) {
            v1 += qPtr[i] * pPtr[i];
            v2 += newQPtr[i] * newPPtr[i];
            v3 += testQ[i] * newPPtr[i];
        }

        cout << "check values: " << v1 << "," << v2 << "," << v3 << endl;

        delete newQ;
        delete newP;
        delete uForQ;
        delete testQ;
    }

}
#endif //VALUEDISTRIBUTION_SVDDISANALYSIS2_H
