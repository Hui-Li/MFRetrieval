#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <iomanip>
#include "Matrix.h"

using namespace std;


double calNorms(const Matrix &m, vector<VectorElement> &norms) {
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

void CalNorms(const vector<int> rowIDs, const Matrix &m, vector<double> &norms) {

    norms.resize(rowIDs.size());

    for (int i = 0; i < rowIDs.size(); i++) {
        int rowID = rowIDs[i];

        double normValue = 0;
        const double *rowPtr = m.getRowPtr(rowID);

        for (int colIndex = 0; colIndex < m.colNum; colIndex++) {
            normValue += rowPtr[colIndex] * rowPtr[colIndex];
        }

        norms[i] = sqrt(normValue);
    }
}

void CalPartitionNorms(const vector<int> rowIDs, const Matrix &m, vector<double> &norms, vector<vector<double> > &partitionNorms, const int partitionNum) {

    norms.resize(rowIDs.size());
    partitionNorms.resize(rowIDs.size(), vector<double>(partitionNum));

    int partitionSize = m.colNum / partitionNum;

    cout << "partitionSize:" << partitionSize << endl;

    for (int i = 0; i < rowIDs.size(); i++) {
        int rowID = rowIDs[i];

        double normValue = 0;
        double partitionNormValue = 0;
        int partitionIndex = 0;
        int innerPartitionIndex = 0;

        const double *rowPtr = m.getRowPtr(rowID);

        for (int colIndex = 0; colIndex < m.colNum; colIndex++) {

            double temp = rowPtr[colIndex] * rowPtr[colIndex];
            normValue += temp;
            partitionNormValue += temp;
            innerPartitionIndex++;

            if (partitionIndex != partitionNum - 1 && innerPartitionIndex == partitionSize) {
                partitionNorms[i][partitionIndex] = sqrt(partitionNormValue);
                partitionIndex++;
                innerPartitionIndex = 0;
                partitionNormValue = 0;
            }
            else if (partitionIndex == partitionNum - 1 && colIndex == m.colNum - 1) {  // Last partition
                partitionNorms[i][partitionIndex] = sqrt(partitionNormValue);
            }
        }

        norms[i] = sqrt(normValue);
    }
}

inline double innerProduct(const double *a, const double *b, const int dim) {
    double value = 0;
    for(int i=0;i<dim;i++){
        value+= a[dim] * b[dim];
    }
    return value;
}

void output(const string outputPath, vector<VectorElement> &norms){

    ofstream outputFile(outputPath.c_str());

    for(int colIndex = 0; colIndex < norms.size(); colIndex++) {
        outputFile << norms[colIndex].data << endl;
    }

//    for(int colIndex = 0; colIndex < norms2.size(); colIndex++) {
//        outputFile << norms2[colIndex].data << endl;
//    }

    outputFile.close();
}

void outputOriginalNorm(string dataset, string type){
    Matrix data;
    data.readData("../../../data/" + dataset + "/" + type + ".txt");
    cout << type << ":" << data.rowNum << "," << data.colNum << endl;

    vector<VectorElement> norms;
    calNorms(data, norms);
    sort(norms.begin(), norms.end(), greater<VectorElement>());
    double max = norms[0].data;
    double min = norms[norms.size()].data;

    double binWidth = 0.2;
    double up = ((int)(max / binWidth)) * binWidth + binWidth;
    int binSize = up/binWidth;

    cout << min << "," << max << "," << up << "," << binSize << "," << binWidth << endl;

    vector<double> bound(binSize+1, 0);
    vector<double> count(binSize, 0);
    for(int i=1;i<binSize+1;i++){
        bound[i] = bound[i-1] + binWidth;
    }

//    for(auto value:bound){
//        cout << value << endl;
//    }
//    exit(0);
    bound[binSize] = max;

    for(int i=0;i<norms.size();i++){
        double norm = norms[i].data;

        for(int j=1;j<bound.size();j++){
            if(norm <= bound[j]){
                count[j-1]++;
                break;
            }
        }
    }

    string outputPath = dataset + "_" + type + "Norms.txt";
    ofstream outputFile(outputPath.c_str());

    for(int i = 0; i < binSize; i++) {
//        outputFile << i + 1 << "," <<(bound[i] + bound[i+1]) / 2 << "," << count[i] << endl;

//        string s = "";
//        stringstream ss;
//        if(i%4==0){
//            ss << fixed << setprecision(1) << bound[i];
//            s = ss.str();
//        }

//        outputFile << s << "," << bound[i] << "," << count[i] << endl;
        outputFile << bound[i] << "," << count[i] << endl;
    }

    outputFile.close();

}

void normRatio(int qSampleSize, int partitionNum){
    Matrix Q;
    Q.readData("../q.txt");
    cout << "Q:" << Q.rowNum << "," << Q.colNum << endl;

    Matrix P;
    P.readData("../p.txt");
    cout << "P:" << P.rowNum << "," << P.colNum << endl;

    random_device rd;     // only used once to initialise (seed) engine
    mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    uniform_int_distribution<int> uni(0,Q.rowNum-1); // guaranteed unbiased

    vector<double> pNorms;
    vector<vector<double> > pPartitionNorms;
    vector<double> qNorms;
    vector<vector<double> > qPartitionNorms;

    vector<int> pIDs(P.rowNum);
    for (int i = 0; i < P.rowNum; i++) {
        pIDs[i] = i;
    }

    vector<int> sampleQIDs(qSampleSize);
    for (int i = 0; i < qSampleSize; i++) {
        sampleQIDs[i] = uni(rng);
    }

    CalPartitionNorms(pIDs, P, pNorms, pPartitionNorms, partitionNum);
    CalPartitionNorms(sampleQIDs, Q, qNorms, qPartitionNorms, partitionNum);

    vector<string> ratio(qSampleSize * P.rowNum);
    int index = 0;
    // calculate ratio
    for (int i = 0; i < sampleQIDs.size(); i++) {
        for (int j = 0; j < P.rowNum; j++) {
            double LB1 = qNorms[i] * pNorms[j];
            double LB2 = 0;

            for (int k = 0; k < partitionNum; k++) {
                LB2 += pPartitionNorms[j][k] * qPartitionNorms[i][k];
            }

            ratio[index] = to_string(LB1/LB2);
            index++;
        }
    }

    string resultPath = "ratio-" + to_string(qSampleSize) + "-" + to_string(partitionNum) + ".txt";
    ofstream file(resultPath.c_str());

    for (auto line:ratio) {
        file << line << endl;
    }
}

void normRatioWithExactDis(int qSampleSize){
    Matrix Q;
    Q.readData("../q.txt");
    cout << "Q:" << Q.rowNum << "," << Q.colNum << endl;

    Matrix P;
    P.readData("../p.txt");
    cout << "P:" << P.rowNum << "," << P.colNum << endl;

    random_device rd;     // only used once to initialise (seed) engine
    mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    uniform_int_distribution<int> uni(0,Q.rowNum-1); // guaranteed unbiased

    vector<double> pNorms;
    vector<double> qNorms;

    vector<int> pIDs(P.rowNum);
    for (int i = 0; i < P.rowNum; i++) {
        pIDs[i] = i;
    }

    vector<int> sampleQIDs(qSampleSize);
    for (int i = 0; i < qSampleSize; i++) {
        sampleQIDs[i] = uni(rng);
    }

    CalNorms(pIDs, P, pNorms);
    CalNorms(sampleQIDs, Q, qNorms);

    vector<string> ratio(qSampleSize * P.rowNum);
    int index = 0;
    // calculate ratio
    for (int i = 0; i < sampleQIDs.size(); i++) {
        for (int j = 0; j < P.rowNum; j++) {
            double LB1 = qNorms[i] * pNorms[j];
            double exact = innerProduct(Q.getRowPtr(sampleQIDs[i]), P.getRowPtr(j), P.colNum);

            ratio[index] = to_string(exact/LB1);
            index++;
        }
    }

    string resultPath = "exact-ratio-" + to_string(qSampleSize) + ".txt";
    ofstream file(resultPath.c_str());

    for (auto line:ratio) {
        file << line << endl;
    }
}

int main() {

    string data[] = {"MovieLens", "Yelp", "Netflix", "KDD"};
    for(int i=0;i<4;i++) {
        outputOriginalNorm(data[i], "q");
        outputOriginalNorm(data[i], "p");
    }
//    int qSampleSize = 100;
//    int partitionNum = 8;
//    normRatio(qSampleSize, partitionNum);
//    normRatioWithExactDis(qSampleSize);

    return 0;
}