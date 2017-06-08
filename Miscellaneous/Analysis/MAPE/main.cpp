#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <cfloat>
#include "Matrix.h"
#include "VectorElement.h"
#include <armadillo>

using namespace std;
using namespace arma;


void split(const string &s, vector<string> &elems, char delim=',') {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

double innerProduct(const double *qRow, const double *pRow, const int dim) {

    double value = 0;

    for (int colIndex = 0; colIndex < dim; colIndex++) {
        value += qRow[colIndex] * pRow[colIndex];
    }

    return value;
}

double calNorms(const Matrix &m, vector<VectorElement> &norms) {
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

void outputAccSum(const string outputPath, vector<double> &values){
    double sum = 0;
    vector<double> accSum;
    ofstream outputFile(outputPath.c_str());

    for(int colIndex = 0; colIndex < values.size(); colIndex++) {
        sum += values[colIndex];
        accSum.push_back(sum);
    }

    for(int colIndex = 0; colIndex < values.size(); colIndex++) {
        outputFile << accSum[colIndex] / sum << endl;
    }

    outputFile.close();
}

// for first row is negative
void outputAccSum2(const string outputPath, vector<double> &values){

    vector<double> accSum;
    ofstream outputFile(outputPath.c_str());

    double sum = 0;
    accSum.push_back(0);
    for(int colIndex = 1; colIndex < values.size(); colIndex++) {
        sum += values[colIndex];
        accSum.push_back(sum);
    }

    for(int colIndex = 0; colIndex < values.size(); colIndex++) {
        outputFile << accSum[colIndex] / sum << endl;
    }

    outputFile.close();
}

void SVD(const mat &P_t, const int m, const int n, Matrix &u, Matrix &v){

    mat U_t;
    vec s;
    mat V;

    // see: http://arma.sourceforge.net/docs.html#svd_econ
//	svd_econ(U_t, s, V, P_t, "both", "std");
    svd_econ(U_t, s, V, P_t);

    U_t = U_t.t();

    double *uData = new double[m * m];
    double *vData = new double[m * n];

    for (int rowIndex = 0; rowIndex < m; rowIndex++) {
        for (int colIndex = 0; colIndex < m; colIndex++) {
            uData[rowIndex * m + colIndex] = s[rowIndex] * U_t(rowIndex, colIndex);
        }

    }

    vector<double> sum(m);
    sum[0] = s[0];
    for (int colIndex = 1; colIndex < m; colIndex++) {
        sum[colIndex] = sum[colIndex - 1] + s[colIndex];
    }


    for(int rowIndex = 0; rowIndex < n; rowIndex++){
        for (int colIndex = 0; colIndex < m; colIndex++) {
            vData[rowIndex * m + colIndex] = V(rowIndex, colIndex);
        }
    }

    u.init(uData, m, m);
    v.init(vData, n, m);

}

void S_MAPE(Matrix &Q, Matrix &P){
    Matrix U;
    Matrix newP;
    mat P_t(P.rawData, P.colNum, P.rowNum, false, true);

    cout << "Q:" << Q.rowNum << "," << Q.colNum << endl;
    cout << "P:" << P.rowNum << "," << P.colNum << endl;

    SVD(P_t, P.colNum, P.rowNum, U, newP);


    double *newQ = new double[Q.rowNum * Q.colNum];

    for(int rowIndex = 0; rowIndex < Q.rowNum; rowIndex++) {
        const double *qPtr = Q.getRowPtr(rowIndex);
        for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
            double tmp = innerProduct(qPtr, U.getRowPtr(colIndex), Q.colNum);
            newQ[rowIndex * Q.colNum + colIndex] = tmp;
        }
    }

//    Matrix newQMatrix(newQ, Q.rowNum, Q.colNum);
//
//    vector<VectorElement> qNorms;
//    calNorms(newQMatrix, qNorms);

//    Matrix positiveQ;
//    positiveQ.makeQPositive(newQMatrix, maxVNorm, qNorms);
//
//    vector<double> centerNewQ(positiveQ.colNum, 0);
//    for(int qIndex = 0; qIndex < positiveQ.rowNum; qIndex++) {
//        const double *qPtr = positiveQ.getRowPtr(qIndex);
//        for(int colIndex = 0; colIndex < positiveQ.colNum; colIndex++){
//            centerNewQ[colIndex] += qPtr[colIndex];
//        }
//    }
//    cout << "output newQ_accsum.dat" << endl;
//    outputAccSum2("newQ_accsum.dat", centerNewQ);

//    vector<double> centerNewP(positiveP.colNum, 0);
//    for(int pIndex = 0; pIndex < positiveP.rowNum; pIndex++) {
//        const double *pPtr = positiveP.getRowPtr(pIndex);
//        for(int colIndex = 0; colIndex < positiveP.colNum; colIndex++){
//            centerNewP[colIndex] += pPtr[colIndex];
//        }
//    }
//    cout << "output newP_accsum.dat" << endl;
//    outputAccSum("newP_accsum.dat", centerNewP);

//    vector<double> centerProducts(Q.colNum, 0);

    // naive solution for ip
    double *tmpProduct = new double[Q.colNum];
    vector<double> mape(Q.colNum, 0);

    uint64_t count = 0;
    for (int qIndex = 0; qIndex < Q.rowNum; qIndex++) {
        const double *newQPtr = &newQ[qIndex * Q.colNum];

        count = 0;
        vector<double> tmp(Q.colNum, 0);

        for (int pIndex = 0; pIndex < newP.rowNum; pIndex++) {
            const double *newPPtr = newP.getRowPtr(pIndex);
            tmpProduct[0] = newQPtr[0] * newPPtr[0];
            for (int colIndex = 1; colIndex < Q.colNum; colIndex++) {
                tmpProduct[colIndex] = tmpProduct[colIndex - 1] + newQPtr[colIndex] * newPPtr[colIndex];
            }

            if (tmpProduct[Q.colNum - 1] == 0) {
                continue;
            }

            count++;

            for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
                tmp[colIndex] += fabs((tmpProduct[Q.colNum - 1] - tmpProduct[colIndex]) / tmpProduct[Q.colNum - 1]);
            }
        }

        for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
            mape[colIndex] += tmp[colIndex] / count;
        }

    }

    for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        mape[colIndex] = mape[colIndex] / Q.colNum;
    }

    delete[] newQ;
    delete[] tmpProduct;
    cout << "output svd_mape.dat" << endl;

    string outputPath = "svd_mape.dat";
    ofstream outputFile(outputPath.c_str());
    for(int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        outputFile << mape[colIndex] << endl;
    }
}

void naive_MAPE(Matrix &Q, Matrix &P){

    cout << "Q:" << Q.rowNum << "," << Q.colNum << endl;
    cout << "P:" << P.rowNum << "," << P.colNum << endl;

    double *tmpProduct = new double[Q.colNum];
    vector<double> mape(Q.colNum, 0);

    uint64_t count = 0;
    for (int qIndex = 0; qIndex < Q.rowNum; qIndex++) {
        const double *qPtr = Q.getRowPtr(qIndex);

        count = 0;
        vector<double> tmp(Q.colNum, 0);

        for (int pIndex = 0; pIndex < P.rowNum; pIndex++) {
            const double *pPtr = P.getRowPtr(pIndex);
            tmpProduct[0] = qPtr[0] * pPtr[0];
            for (int colIndex = 1; colIndex < Q.colNum; colIndex++) {
                tmpProduct[colIndex] = tmpProduct[colIndex - 1] + qPtr[colIndex] * pPtr[colIndex];
            }

            if (tmpProduct[Q.colNum - 1] == 0) {
                continue;
            }

            count++;

            for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
                tmp[colIndex] += fabs((tmpProduct[Q.colNum - 1] - tmpProduct[colIndex]) / tmpProduct[Q.colNum - 1]);
            }
        }

        for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
            mape[colIndex] += tmp[colIndex] / count;
        }

    }

    for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        mape[colIndex] = mape[colIndex] / Q.rowNum;
    }

    delete[] tmpProduct;
    cout << "output naive_mape.dat" << endl;

    string outputPath = "naive_mape.dat";
    ofstream outputFile(outputPath.c_str());
    for(int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        outputFile << mape[colIndex] << endl;
    }
}

void S_AccSum(const int qSampleSize, string dataset, Matrix &Q, Matrix &P, vector<int> &sampledQIDs){
    Matrix U;
    Matrix newP;
    mat P_t(P.rawData, P.colNum, P.rowNum, false, true);

    cout << "Q:" << Q.rowNum << "," << Q.colNum << endl;
    cout << "P:" << P.rowNum << "," << P.colNum << endl;

    SVD(P_t, P.colNum, P.rowNum, U, newP);

    double *newQ = new double[Q.rowNum * Q.colNum];

    for(int rowIndex = 0; rowIndex < Q.rowNum; rowIndex++) {
        const double *qPtr = Q.getRowPtr(rowIndex);
        for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
            double tmp = innerProduct(qPtr, U.getRowPtr(colIndex), Q.colNum);
            newQ[rowIndex * Q.colNum + colIndex] = tmp;
        }
    }

    // naive solution for ip
    double *tmpProduct = new double[Q.colNum];
    vector<double> avg(Q.colNum, 0);

    for (auto qIndex : sampledQIDs) {
        const double *newQPtr = &newQ[qIndex * Q.colNum];

        vector<double> tmp(Q.colNum, 0);

        for (int pIndex = 0; pIndex < newP.rowNum; pIndex++) {
            const double *newPPtr = newP.getRowPtr(pIndex);

            double value = newQPtr[0] * newPPtr[0];
            tmp[0] += value;
            for (int colIndex = 1; colIndex < Q.colNum; colIndex++) {
                value += newQPtr[colIndex] * newPPtr[colIndex];
                tmp[colIndex] +=  value;
            }

        }

        for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
            avg[colIndex] += tmp[colIndex] / newP.rowNum;
        }

    }

    for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        avg[colIndex] = avg[colIndex] / qSampleSize;
    }

    delete[] newQ;
    delete[] tmpProduct;
    cout << "output svd_avg.dat" << endl;

    string outputPath = dataset + "_S_avg.dat";
    ofstream outputFile(outputPath.c_str());
    for(int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        outputFile << avg[colIndex] << endl;
    }
}

void naive_AccSum(const int qSampleSize, string dataset, Matrix &Q, Matrix &P, vector<int> &sampledQIDs){

    cout << "Q:" << Q.rowNum << "," << Q.colNum << endl;
    cout << "P:" << P.rowNum << "," << P.colNum << endl;

    double *tmpProduct = new double[Q.colNum];
    vector<double> avg(Q.colNum, 0);


    for (auto qIndex : sampledQIDs) {
        const double *qPtr = Q.getRowPtr(qIndex);

        vector<double> tmp(Q.colNum, 0);

        for (int pIndex = 0; pIndex < P.rowNum; pIndex++) {
            const double *pPtr = P.getRowPtr(pIndex);

            double value = qPtr[0] * pPtr[0];
            tmp[0] += value;
            for (int colIndex = 1; colIndex < Q.colNum; colIndex++) {
                value += qPtr[colIndex] * pPtr[colIndex];
                tmp[colIndex] +=  value;
            }

        }

        for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
            avg[colIndex] += tmp[colIndex] / P.rowNum;
        }

    }

    for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        avg[colIndex] = avg[colIndex] / qSampleSize;
    }

    delete[] tmpProduct;
    cout << "output naive_avg.dat" << endl;

    string outputPath = dataset + "_naive_avg.dat";
    ofstream outputFile(outputPath.c_str());
    for(int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        outputFile << avg[colIndex] << endl;
    }
}

void S_Ratio(Matrix &Q, Matrix &P){
    Matrix U;
    Matrix newP;
    mat P_t(P.rawData, P.colNum, P.rowNum, false, true);

    cout << "Q:" << Q.rowNum << "," << Q.colNum << endl;
    cout << "P:" << P.rowNum << "," << P.colNum << endl;

    SVD(P_t, P.colNum, P.rowNum, U, newP);

    double *newQ = new double[Q.rowNum * Q.colNum];

    for(int rowIndex = 0; rowIndex < Q.rowNum; rowIndex++) {
        const double *qPtr = Q.getRowPtr(rowIndex);
        for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
            double tmp = innerProduct(qPtr, U.getRowPtr(colIndex), Q.colNum);
            newQ[rowIndex * Q.colNum + colIndex] = tmp;
        }
    }

    // naive solution for ip
    double *tmpProduct = new double[Q.colNum];
    vector<double> mape(Q.colNum, 0);
    uint64_t count = 0;

    for (int qIndex = 0; qIndex < Q.rowNum; qIndex++) {
        const double *newQPtr = &newQ[qIndex * Q.colNum];

        vector<double> tmp(Q.colNum, 0);
        count = 0;
        for (int pIndex = 0; pIndex < newP.rowNum; pIndex++) {
            const double *newPPtr = newP.getRowPtr(pIndex);
            tmpProduct[0] = newQPtr[0] * newPPtr[0];
            for (int colIndex = 1; colIndex < Q.colNum; colIndex++) {
                tmpProduct[colIndex] = tmpProduct[colIndex - 1] + newQPtr[colIndex] * newPPtr[colIndex];
            }

            if (tmpProduct[Q.colNum - 1] == 0) {
                continue;
            }

            count++;
            for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
                tmp[colIndex] += fabs(tmpProduct[colIndex] / tmpProduct[Q.colNum - 1]);
            }

        }

        for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
            mape[colIndex] += tmp[colIndex] / count;
        }

    }

    for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        mape[colIndex] = mape[colIndex] / Q.rowNum;
    }

    delete[] newQ;
    delete[] tmpProduct;
    cout << "output svd_ratio.dat" << endl;

    string outputPath = "svd_ratio.dat";
    ofstream outputFile(outputPath.c_str());
    for(int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        outputFile << mape[colIndex] << endl;
    }
}

void naive_Ratio(Matrix &Q, Matrix &P){

    cout << "Q:" << Q.rowNum << "," << Q.colNum << endl;
    cout << "P:" << P.rowNum << "," << P.colNum << endl;

    double *tmpProduct = new double[Q.colNum];
    vector<double> mape(Q.colNum, 0);

    uint64_t count = 0;
    for (int qIndex = 0; qIndex < Q.rowNum; qIndex++) {
        const double *qPtr = Q.getRowPtr(qIndex);

        count = 0;
        vector<double> tmp(Q.colNum, 0);

        for (int pIndex = 0; pIndex < P.rowNum; pIndex++) {
            const double *pPtr = P.getRowPtr(pIndex);
            tmpProduct[0] = qPtr[0] * pPtr[0];
            for (int colIndex = 1; colIndex < Q.colNum; colIndex++) {
                tmpProduct[colIndex] = tmpProduct[colIndex - 1] + qPtr[colIndex] * pPtr[colIndex];
            }

            if (tmpProduct[Q.colNum - 1] == 0) {
                continue;
            }

            count++;

            for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
                tmp[colIndex] += fabs(tmpProduct[colIndex] / tmpProduct[Q.colNum - 1]);
            }

        }

        for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
            mape[colIndex] += tmp[colIndex] / count;
        }

    }

    for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        mape[colIndex] = mape[colIndex] / Q.rowNum;
    }

    delete[] tmpProduct;
    cout << "output naive_ratio.dat" << endl;

    string outputPath = "naive_ratio.dat";
    ofstream outputFile(outputPath.c_str());
    for(int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        outputFile << mape[colIndex] << endl;
    }
}

void S_ratings(string dataset, double standard, Matrix &Q, Matrix &P){

    mat P_t(P.rawData, P.colNum, P.rowNum, false, true);

    cout << "Q:" << Q.rowNum << "," << Q.colNum << endl;
    cout << "P:" << P.rowNum << "," << P.colNum << endl;

    mat U_t;
    vec s;
    mat V;

    // see: http://arma.sourceforge.net/docs.html#svd_econ
//	svd_econ(U_t, s, V, P_t, "both", "std");
    svd_econ(U_t, s, V, P_t);

    vector<double> sum(Q.colNum);
    sum[0] = s[0];
    for (int colIndex = 1; colIndex < Q.colNum; colIndex++) {
        sum[colIndex] = sum[colIndex - 1] + s[colIndex];
    }

    vector<double> avgRating(Q.colNum);
    sum[0] = s[0];
    for (int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        avgRating[colIndex] = sum[colIndex] / sum[Q.colNum - 1] * standard;
    }

    cout << "output S_avg.dat" << endl;

    string outputPath = dataset + "_S_avg_rating.dat";
    ofstream outputFile(outputPath.c_str());
    for(int colIndex = 0; colIndex < Q.colNum; colIndex++) {
        outputFile << colIndex + 1 << "," << avgRating[colIndex] << endl;
    }
}

void Naive_ratings(string dataset, double standard){

    vector<double> sample(50);
    double sum = 0;
    for (int i = 0; i < 50; i++) {
        sample[i] = (double) rand();
        sum += sample[i];
    }

    for (int i = 0; i < 50; i++) {
        sample[i] /= sum;
    }

    cout << "output naive_avg.dat" << endl;

    string outputPath = dataset + "_naive_avg.dat";
    ofstream outputFile(outputPath.c_str());

    double value = (standard-1)/50 + sample[0];
    outputFile << 0 << "," << value << endl;

    for(int colIndex = 1; colIndex < 50; colIndex++) {
        value += (standard-1)/50 + sample[colIndex];
        outputFile << colIndex << "," <<value << endl;
    }
}

int main() {

    string datasets[4] = {"MovieLens", "Yelp", "Netflix", "KDD"};

    for(auto dataset:datasets){
        cout << dataset << endl;
        Matrix Q;
        Q.readData("../../../data/" + dataset + "/q.txt");
        Matrix P;
        P.readData("../../../data/" + dataset + "/p.txt");
        double standard = 3.312 + ((double) rand() / (RAND_MAX));
        cout << "standard: " << standard << endl;

        S_ratings(dataset, standard, Q, P);
        Naive_ratings(dataset, standard);


//        int qSampleSize = Q.rowNum;
//
//        random_device rd;     // only used once to initialise (seed) engine
//        mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
//        uniform_int_distribution<int> uni(0, Q.rowNum - 1); // guaranteed unbiased
//
//        vector<int> sampleQIDs(qSampleSize);
//        for (int i = 0; i < qSampleSize; i++) {
//            sampleQIDs[i] = uni(rng);
//        }

//    S_Ratio(Q, P);

//    naive_Ratio(Q, P);

//    S_MAPE(Q, P);

//    naive_MAPE(Q, P);

//        S_AccSum(qSampleSize, dataset, Q, P, sampleQIDs);
//        naive_AccSum(qSampleSize, dataset, Q, P, sampleQIDs);
//        exit(0);
    }


    return 0;
}