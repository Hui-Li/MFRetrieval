#include <iostream>
#include "structs/Matrix.h"
#include "util/Monitor.h"
#include "lsh/psdLSH.h"

int main() {

    string qFilePath = "../data/q.txt";
    string pFilePath = "../data/p.txt";
    string resultPath = "result.txt";

    Monitor timer;

    Matrix Q, P, transformedQ, transformedP;
    Q.readData(qFilePath);
    P.readData(pFilePath);

    timer.start();
    transformedQ.transformQMatrix(Q);
    transformedP.transformPMatrix(P);
    timer.stop();

    cout << "time for data preparing: " << timer.getElapsedTime() << " secs" << endl;

    cout << "Q: " << Q.rowNum << "," << Q.colNum << endl;
    cout << "P: " << P.rowNum << "," << P.colNum << endl;

    timer.start();

    /////////////////////// Parameters ////////////////////////////
    int K = 1;
    psdLSH::Parameter param;
    param.hash_table_size = 521;
    param.num_of_hash_tables = 30;
    param.D = transformedP.colNum;
    param.T = GAUSSIAN;
    param.windows_size = 10;
    /////////////////////// Parameters ////////////////////////////

    psdLSH lsh;
    lsh.init(param);
    lsh.hash(transformedP);

    timer.stop();

    cout << "time for indexing: " << timer.getElapsedTime() << " secs" << endl;

    Matrix::Accessor accessor(transformedP);
    Metric metric(transformedP.colNum, L2_DIST);

    vector<Scanner<Matrix::Accessor> > scanners(transformedQ.rowNum, Scanner<Matrix::Accessor>(accessor, metric, K));

    timer.start();

    for (int i = 0; i != transformedQ.rowNum; ++i) {
        Scanner<Matrix::Accessor> &scanner = scanners[i];
        lsh.query(transformedQ[i], scanner);
    }

    timer.stop();

    cout << "time for query: " << timer.getElapsedTime() << " secs" << endl;

    ofstream file(resultPath.c_str());

    for (int qIndex = 0; qIndex < transformedQ.rowNum; qIndex++) {
        Scanner<Matrix::Accessor> &scanner = scanners[qIndex];
        std::vector<std::pair<double, unsigned> > &results = scanner.topk().getTopk();
        for (int j = 0; j < K; j++) {
            file << qIndex << "," << results[j].second << "," << results[j].first << endl;
        }
    }
    file.close();

    return 0;
}