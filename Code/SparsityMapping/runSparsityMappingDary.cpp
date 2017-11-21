#include "util/Base.h"
#include "structs/Matrix.h"
#include "util/Monitor.h"
#include "alg/TernaryBaseSchema.h"
#include "util/FileUtil.h"
#include "alg/DaryBaseSchema.h"

int main() {

    //////////////// Parameters /////////////////
    int k = 10;
    int d = 100;
    // input must be normalized vectors
    string q_file = "../Test_Data/q.txt";
    string p_file = "../Test_Data/p.txt";
    string output_path = "./results/Dary-results.txt";
    //////////////// Parameters /////////////////

    Matrix Q, P;
    Q.readData(q_file);
    P.readData(p_file);

    cout << "Q: row " << Q.rowNum << ", col " << Q.colNum << endl;
    cout << "P: row " << P.rowNum << ", col " << P.colNum << endl;

    Monitor timer;

    timer.start();

    DaryBaseSchema dbSchema(Q, P, d);
    dbSchema.data_prepare();

    timer.stop();

    cout << "time for data preparing: " << timer.getElapsedTime() << " secs" << endl;

    timer.start();
    VectorElement *results = new VectorElement[Q.rowNum * k];
    dbSchema.query_topk(k, results);

    timer.stop();

    cout << "time for query: " << timer.getElapsedTime() << " secs" << endl;

    FileUtil::outputResult(k, P, Q, results, output_path);

    delete[] results;

    return 0;
}