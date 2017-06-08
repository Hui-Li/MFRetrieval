#ifndef VALUEDISTRIBUTION_SVDDISANALYSIS_H
#define VALUEDISTRIBUTION_SVDDISANALYSIS_H

#include <armadillo>
#include "Matrix.h"
using namespace std;
using namespace arma;

namespace SVDDisAnalysis {

    void SVD(const mat &P_t, const int m, const int n, Matrix &Q, Matrix &newQ, Matrix &newP, Matrix &u) {

        mat U_t;
        vec s;
        mat V;

        svd_econ(U_t, s, V, P_t);

        U_t = U_t.t();

        double *uData = new double[m * m];
        double *vData = new double[m * n];

        for (int rowIndex = 0; rowIndex < m; rowIndex++) {
            for (int colIndex = 0; colIndex < m; colIndex++) {
                uData[rowIndex * m + colIndex] = s[rowIndex] * U_t(rowIndex, colIndex);
            }
        }

        u.init(uData, m, m);

        for (int QIndex = 0; QIndex < newQ.rowNum; QIndex++) {

            double *newQPtr = newQ.getRowPtr(QIndex);
            double *qPtr = Q.getRowPtr(QIndex);
            for (int rowIndex = u.rowNum - 1; rowIndex >= 0; rowIndex--) {
                newQPtr[rowIndex] = 0;
                const double *uPtr = u.getRowPtr(rowIndex);

                for (int colIndex = 0; colIndex < u.colNum; colIndex++) {
                    newQPtr[rowIndex] += uPtr[colIndex] * qPtr[colIndex];
                }
            }
        }

        for(int rowIndex = 0; rowIndex < n; rowIndex++){
            for (int colIndex = 0; colIndex < m; colIndex++) {
                vData[rowIndex * m + colIndex] = V(rowIndex, colIndex);
            }
        }

        newP.init(vData, n, m);

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
        Matrix *u = new Matrix();
        SVD(P_t, P->colNum, P->rowNum, *Q, *newQ, *newP, *u);

        vector<double> avg(newP->colNum, 0);
        for (int rowIndex = 0; rowIndex < newP->rowNum; rowIndex++) {
            const double *ptr = newP->getRowPtr(rowIndex);
            for (int colIndex = 0; colIndex < newP->colNum; colIndex++) {
                avg[colIndex] += fabs(ptr[colIndex]);
            }
        }

        for (int colIndex = 0; colIndex < newP->colNum; colIndex++) {
            avg[colIndex] /= newP->rowNum;
        }

        cout << "P:" << newP->rowNum << "," << newP->colNum << endl;
        string pFile = data + "_pSVDAvg.txt";
        ofstream outputFile(pFile.c_str());

        for(int colIndex = 0; colIndex < newP->colNum; colIndex++) {
            outputFile << colIndex + 1 << "," << avg[colIndex] << endl;
        }

        outputFile.close();

        vector<double> avg2(newQ->colNum, 0);
        for (int rowIndex = 0; rowIndex < newQ->rowNum; rowIndex++) {
            const double *qtr = newQ->getRowPtr(rowIndex);
            for (int colIndex = 0; colIndex < newQ->colNum; colIndex++) {
                avg2[colIndex] += fabs(qtr[colIndex]);
            }
        }

        for (int colIndex = 0; colIndex < newQ->colNum; colIndex++) {
            avg2[colIndex] /= newQ->rowNum;
        }

        cout << "Q:" << newQ->rowNum << "," << newQ->colNum << endl;
        string qFile = data + "_qSVDAvg.txt";
        ofstream outputFile2(qFile.c_str());

        for(int colIndex = 0; colIndex < newQ->colNum; colIndex++) {
            outputFile2 << colIndex + 1 << "," << avg2[colIndex] << endl;
        }

        outputFile2.close();

        const double *qPtr = Q->getRowPtr(131);
        const double *pPtr = P->getRowPtr(444);

        const double *newQPtr = newQ->getRowPtr(131);
        const double *newPPtr = newP->getRowPtr(444);

        double *testQ = new double[Q->colNum];
        transferQ(u, qPtr, testQ);

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
        delete u;
        delete testQ;
    }

    void calNorms(const Matrix &m, vector<VectorElement> &norms) {
        norms.resize(m.rowNum);

        for (int rowID = 0; rowID < m.rowNum; rowID++) {
            double norm = 0;
            const double *row = m.getRowPtr(rowID);
            for (int colIndex = 0; colIndex < m.colNum; colIndex++) {
                norm += row[colIndex] * row[colIndex];
            }
            norm = sqrt(norm);
            norms[rowID] = VectorElement(rowID, norm);
        }

    }

    double getSVDQBinWidth(string data){
        if(data=="MovieLens"){
            return 50;
        } else if (data=="Yelp"){
            return 100;
        } else if (data=="Netflix"){
            return 50;
        } else if (data=="KDD"){
            return 200;
        }
        return -1;
    }

    double getSVDPBinWidth(string data){
        if(data=="MovieLens"){
            return 0.04;
        } else if (data=="Yelp"){
            return 0.01;
        } else if (data=="Netflix"){
            return 0.01;
        } else if (data=="KDD"){
            return 0.002;
        }
        return -1;
    }

    void analysisForLength(string data){
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
        Matrix *u = new Matrix();
        SVD(P_t, P->colNum, P->rowNum, *Q, *newQ, *newP, *u);

        // calculate length for each
        vector<VectorElement> qNorms;
        calNorms(*newQ, qNorms);
        sort(qNorms.begin(), qNorms.end(), greater<VectorElement>());

        vector<VectorElement> pNorms;
        calNorms(*newP, pNorms);
        sort(pNorms.begin(), pNorms.end(), greater<VectorElement>());

        double qMax = qNorms[0].data;
        double qMin = qNorms[qNorms.size()].data;

        double qBinWidth = getSVDQBinWidth(data);
        double qUP = ((int)(qMax / qBinWidth)) * qBinWidth + qBinWidth;
        int qBinSize = qUP/qBinWidth;

        cout << qMin << "," << qMax << "," << qUP << "," << qBinSize << "," << qBinWidth << endl;

        vector<double> qBound(qBinSize + 1, 0);
        vector<double> qCount(qBinSize, 0);
        for (int i = 1; i < qBinSize + 1; i++) {
            qBound[i] = qBound[i - 1] + qBinWidth;
        }

        qBound[qBinSize] = qMax;

        for(int i=0;i<qNorms.size();i++){
            double norm = qNorms[i].data;

            for(int j=1;j<qBound.size();j++){
                if(norm <= qBound[j]){
                    qCount[j-1]++;
                    break;
                }
            }
        }

        string outputPath = data + "_qSVDNorms.txt";
        ofstream outputFile(outputPath.c_str());

        for(int i = 0; i < qBinSize; i++) {
            outputFile << qBound[i] << "," << qCount[i] << endl;
        }

        outputFile.close();

        double pMax = pNorms[0].data;
        double pMin = pNorms[pNorms.size()].data;

        double pBinWidth = getSVDPBinWidth(data);
        double pUP = ((int)(pMax / pBinWidth)) * pBinWidth + pBinWidth;
        int pBinSize = pUP/pBinWidth;

        cout << pMin << "," << pMax << "," << pUP << "," << pBinSize << "," << pBinWidth << endl;

        vector<double> pBound(pBinSize + 1, 0);
        vector<double> pCount(pBinSize, 0);
        for (int i = 1; i < pBinSize + 1; i++) {
            pBound[i] = pBound[i - 1] + pBinWidth;
        }

        pBound[pBinSize] = pMax;

        for (int i = 0; i < pNorms.size(); i++) {
            double norm = pNorms[i].data;

            for (int j = 1; j < pBound.size(); j++) {
                if (norm <= pBound[j]) {
                    pCount[j - 1]++;
                    break;
                }
            }
        }

        string outputPath2 = data + "_pSVDNorms.txt";
        ofstream outputFile2(outputPath2.c_str());

        for(int i = 0; i < pBinSize; i++) {
            outputFile2 << pBound[i] << "," << pCount[i] << endl;
        }

        outputFile2.close();

        delete Q;
        delete P;
        delete newQ;
        delete newP;
        delete u;


    }

}
#endif //VALUEDISTRIBUTION_SVDDISANALYSIS_H
