#include "util/Base.h"
#include "util/FileUtil.h"
#include "util/Monitor.h"
#include "util/Calculator.h"
#include "structs/FastHeap.h"

/**
 * Naive solution.
 */
void naive(const int k, const Matrix &q, const Matrix &p, VectorElement *results) {

    Monitor t;
    t.start();

    for (int qID = 0; qID < q.rowNum; qID++) {

        if(qID % 1000 == 0){
            cout << qID << " of " << q.rowNum << endl;
        }

        int heapCount = 0;
        const double *qRow = q[qID];
        VectorElement *heap = results + qID * k;

        for (int pID = 0; pID < k; pID++) {
            const double *pRow = p[pID];
            double value = inner_product(qRow, pRow, q.colNum);
            heap_enqueue(value, pID, heap, &heapCount);
        }

        for (int pID = k; pID < p.rowNum; pID++) {
            const double *pRow = p[pID];
            double value = inner_product(qRow, pRow, q.colNum);

            if (value > heap[0].data) {
                heap_dequeue(heap, &heapCount);
                heap_enqueue(value, pID, heap, &heapCount);
            }
        }
    }

    t.stop();

    cout << "time for query: " << t.getElapsedTime() << " secs" << endl;
}


int main() {

    //////////////// Parameters /////////////////
    int k = 10;
    string q_file = "../data/Test/q.txt";
    string p_file = "../data/Test/p.txt";
    string output_path = "./results/Naive-results.txt";
    //////////////// Parameters /////////////////

    Matrix Q, P;
    Q.readData(q_file);
    P.readData(p_file);

    cout << "Q: row " << Q.rowNum << ", col " << Q.colNum << endl;
    cout << "P: row " << P.rowNum << ", col " << P.colNum << endl;

    VectorElement *results = new VectorElement[Q.rowNum * k];

    naive(k, Q, P, results);

    FileUtil::outputNaiveResult(k, P, Q, results, output_path);

    delete[] results;
}