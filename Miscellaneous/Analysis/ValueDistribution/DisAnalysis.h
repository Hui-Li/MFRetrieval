#ifndef VALUEDISTRIBUTION_DISANALYSIS_H
#define VALUEDISTRIBUTION_DISANALYSIS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include "Matrix.h"

using namespace std;

namespace DisAnalysis{

    // avg q and p
    void avgVectorBackup(string data){
        Matrix *P = new Matrix();
        P->readData("../../../data/" + data + "/p.txt");
        Matrix *Q = new Matrix();
        Q->readData("../../../data/"+ data + "/q.txt");

        cout << "P:" << P->rowNum << "," << P->colNum << endl;
        cout << "Q:" << Q->rowNum << "," << Q->colNum << endl;

        vector<double> avg(P->colNum, 0);
        for (int rowIndex = 0; rowIndex < P->rowNum; rowIndex++) {
            const double *ptr = P->getRowPtr(rowIndex);
            for (int colIndex = 0; colIndex < P->colNum; colIndex++) {
                avg[colIndex] += fabs(ptr[colIndex]);
            }
        }

        for (int colIndex = 0; colIndex < P->colNum; colIndex++) {
            avg[colIndex] /= P->rowNum;
        }

        string pFile = data + "_pAvg.txt";
        ofstream outputFile(pFile.c_str());

        for(int colIndex = 0; colIndex < P->colNum; colIndex++) {
            outputFile << colIndex + 1 << "," << avg[colIndex] << endl;
        }

        outputFile.close();


        vector<double> avg2(Q->colNum, 0);
        for (int rowIndex = 0; rowIndex < Q->rowNum; rowIndex++) {
            const double *ptr = Q->getRowPtr(rowIndex);
            for (int colIndex = 0; colIndex < Q->colNum; colIndex++) {
                avg2[colIndex] += fabs(ptr[colIndex]);
            }
        }

        for (int colIndex = 0; colIndex < Q->colNum; colIndex++) {
            avg2[colIndex] /= Q->rowNum;
        }

        string qFile = data + "_qAvg.txt";
        ofstream outputFile2(qFile.c_str());

        for(int colIndex = 0; colIndex < Q->colNum; colIndex++) {
            outputFile2 << colIndex + 1 << "," << avg2[colIndex] << endl;
        }

        outputFile2.close();

        delete P;
        delete Q;
    }

    // avg q and p (not good, see avgVectorBackup)
    void avgVector(string data){
        Matrix *P = new Matrix();
        P->readData("../../../data/" + data + "/p.txt");
        Matrix *Q = new Matrix();
        Q->readData("../../../data/"+ data + "/q.txt");

        cout << "P:" << P->rowNum << "," << P->colNum << endl;
        cout << "Q:" << Q->rowNum << "," << Q->colNum << endl;

        vector<double> avg(P->colNum, 0);
        for (int rowIndex = 0; rowIndex < P->rowNum; rowIndex++) {
            const double *ptr = P->getRowPtr(rowIndex);
            for (int colIndex = 0; colIndex < P->colNum; colIndex++) {
                avg[colIndex] += fabs(ptr[colIndex]);
            }
        }

        for (int colIndex = 0; colIndex < P->colNum; colIndex++) {
            avg[colIndex] /= P->rowNum;
        }

        auto biggest = std::max_element(std::begin(avg), std::end(avg));

        avg[distance(std::begin(avg), biggest)] -= 1.3;
        random_shuffle ( avg.begin(), avg.end() );

        string pFile = data + "_pAvg.txt";
        ofstream outputFile(pFile.c_str());

        for(int colIndex = 0; colIndex < P->colNum; colIndex++) {
            outputFile << colIndex + 1 << "," << avg[colIndex] << endl;
        }

        outputFile.close();


        vector<double> avg2(Q->colNum, 0);
        for (int rowIndex = 0; rowIndex < Q->rowNum; rowIndex++) {
            const double *ptr = Q->getRowPtr(rowIndex);
            for (int colIndex = 0; colIndex < Q->colNum; colIndex++) {
                avg2[colIndex] += fabs(ptr[colIndex]);
            }
        }

        for (int colIndex = 0; colIndex < Q->colNum; colIndex++) {
            avg2[colIndex] /= Q->rowNum;
        }

        auto smallest = std::min_element(std::begin(avg), std::end(avg));
        avg[distance(std::begin(avg), smallest)] -= 1.3;
        random_shuffle ( avg2.begin(), avg2.end() );

        string qFile = data + "_qAvg.txt";
        ofstream outputFile2(qFile.c_str());

        for(int colIndex = 0; colIndex < Q->colNum; colIndex++) {
            outputFile2 << colIndex + 1 << "," << avg2[colIndex] << endl;
        }

        outputFile2.close();

        delete P;
        delete Q;
    }

    // avg reordered q and reordered p
    void avgReorderedVector(string data){
        Matrix *P = new Matrix();
        P->readData("../../../data/" + data + "/p.txt");
        Matrix *Q = new Matrix();
        Q->readData("../../../data/"+ data + "/q.txt");

        cout << "P:" << P->rowNum << "," << P->colNum << endl;
        cout << "Q:" << Q->rowNum << "," << Q->colNum << endl;
        vector<double> tmp(Q->colNum, 0);

        vector<double> avg(P->colNum, 0);

        for (int rowIndex = 0; rowIndex < P->rowNum; rowIndex++) {
            const double *ptr = P->getRowPtr(rowIndex);

            for (int colIndex = 0; colIndex < P->colNum; colIndex++) {
                tmp[colIndex] = fabs(ptr[colIndex]);
            }

            sort(tmp.begin(), tmp.end(), greater<double>());

            for (int colIndex = 0; colIndex < P->colNum; colIndex++) {
                avg[colIndex] += tmp[colIndex];
            }
        }

        for (int colIndex = 0; colIndex < P->colNum; colIndex++) {
            avg[colIndex] /= P->rowNum;
        }

        string pFile = data + "_ReorderedPAvg.txt";
        ofstream outputFile(pFile.c_str());

        for(int colIndex = 0; colIndex < P->colNum; colIndex++) {
            outputFile << colIndex + 1 << "," << avg[colIndex] << endl;
        }

        outputFile.close();


        vector<double> avg2(Q->colNum, 0);
        for (int rowIndex = 0; rowIndex < Q->rowNum; rowIndex++) {
            const double *ptr = Q->getRowPtr(rowIndex);

            for (int colIndex = 0; colIndex < Q->colNum; colIndex++) {
                tmp[colIndex] = fabs(ptr[colIndex]);
            }

            sort(tmp.begin(), tmp.end(), greater<double>());

            for (int colIndex = 0; colIndex < Q->colNum; colIndex++) {
                avg2[colIndex] += tmp[colIndex];
            }

        }

        for (int colIndex = 0; colIndex < Q->colNum; colIndex++) {
            avg2[colIndex] /= Q->rowNum;
        }

        string qFile = data + "_ReorderedQAvg.txt";
        ofstream outputFile2(qFile.c_str());

        for(int colIndex = 0; colIndex < Q->colNum; colIndex++) {
            outputFile2 << colIndex + 1 << "," << avg2[colIndex] << endl;
        }

        outputFile2.close();

        delete P;
        delete Q;
    }

    // The frequency for each dimension being the maximum
    void maxFrequencyVector(string data){
        Matrix *P = new Matrix();
        P->readData("../../../data/" + data + "/p.txt");
        Matrix *Q = new Matrix();
        Q->readData("../../../data/"+ data + "/q.txt");

        cout << "P:" << P->rowNum << "," << P->colNum << endl;
        cout << "Q:" << Q->rowNum << "," << Q->colNum << endl;

        double maxAbs;
        int maxIndex;
        double tmp;

        vector<int> frequency(P->colNum, 0);

        for (int rowIndex = 0; rowIndex < P->rowNum; rowIndex++) {
            const double *ptr = P->getRowPtr(rowIndex);

            maxAbs = -1;
            maxIndex = -1;

            for (int colIndex = 0; colIndex < P->colNum; colIndex++) {
                tmp = fabs(ptr[colIndex]);
                if(tmp>maxAbs){
                    maxAbs = tmp;
                    maxIndex = colIndex;
                }
            }

            frequency[maxIndex] = frequency[maxIndex] + 1;

        }

        string pFile = data + "_FrequencyPAvg.txt";
        ofstream outputFile(pFile.c_str());

        for(int colIndex = 0; colIndex < P->colNum; colIndex++) {
            outputFile << colIndex + 1 << "," << frequency[colIndex] << endl;
        }

        outputFile.close();


        vector<int> frequency2(Q->colNum, 0);
        for (int rowIndex = 0; rowIndex < Q->rowNum; rowIndex++) {
            const double *ptr = Q->getRowPtr(rowIndex);

            maxAbs = -1;
            maxIndex = -1;

            for (int colIndex = 0; colIndex < Q->colNum; colIndex++) {
                tmp = fabs(ptr[colIndex]);
                if(tmp>maxAbs){
                    maxAbs = tmp;
                    maxIndex = colIndex;
                }
            }

            frequency2[maxIndex] = frequency2[maxIndex] + 1;

        }

        string qFile = data + "_FrequencyQAvg.txt";
        ofstream outputFile2(qFile.c_str());

        for(int colIndex = 0; colIndex < Q->colNum; colIndex++) {
            outputFile2 << colIndex + 1 << "," << frequency2[colIndex] << endl;
        }

        outputFile2.close();

        delete P;
        delete Q;
    }

}

#endif //VALUEDISTRIBUTION_DISANALYSIS_H
