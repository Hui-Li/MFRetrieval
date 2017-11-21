#include "util/Base.h"
#include "structs/Matrix.h"
#include "util/Monitor.h"
#include "alg/TernaryBaseSchema.h"
#include "util/FileUtil.h"

int main() {

    //////////////// Parameters /////////////////
    int k = 10;
    string q_file = "../data/Test/q.txt";
    string p_file = "../data/Test/p.txt";
    string output_path = "./results/Ternary-results.txt";
    //////////////// Parameters /////////////////

    Matrix Q, P;
    Q.readData(q_file);
    P.readData(p_file);

    cout << "Q: row " << Q.rowNum << ", col " << Q.colNum << endl;
    cout << "P: row " << P.rowNum << ", col " << P.colNum << endl;

    int p = Q.colNum * 3; // padding dimensions

    Monitor timer;

    timer.start();

    TernaryBaseSchema tbSchema(Q, P, p);
    tbSchema.data_prepare();

    timer.stop();

    cout << "time for data preparing: " << timer.getElapsedTime() << " secs" << endl;

    timer.start();
    VectorElement *results = new VectorElement[Q.rowNum * k];
    tbSchema.query_topk(k, results);

    timer.stop();

    cout << "time for query: " << timer.getElapsedTime() << " secs" << endl;

    FileUtil::outputResult(k, P, Q, results, output_path);

    delete[] results;

    return 0;
}