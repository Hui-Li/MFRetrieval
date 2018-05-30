#ifndef WEDGENNSAMPLING_H
#define WEDGENNSAMPLING_H

#include <random>
#include "../structs/VectorElement.h"
#include "../structs/Node.h"
#include "../structs/Edge.h"
#include "../util/Base.h"
#include "../util/FileUtil.h"
#include "../util/Calculator.h"
#include "../util/AliasSamplingGenerator.h"
#include "../structs/FastHeap.h"
#include "../util/Monitor.h"
#include "../util/RandomUtil.h"
#include "../util/EvalUtil.h"


class WedgeNNSampling{

private:

    int QNum, PNum, d;
    double *QData = nullptr;
    double *PData = nullptr;

    /**
     * Compute cumulative sums and column norms of dense matrix
     * @param A
     * @param A_row_num
     * @param A_col_num
     * @param colnorms
     * @param cumsums
     */
    inline void cumulative_sums(const double *A, int A_row_num, int A_col_num, double *colnorms, double *cumsums) {
        cumsums[0] = 0.0;
        for (int col_index = 0; col_index < A_col_num; col_index++) {
            colnorms[col_index] = 0.0;
            for (int row_index = col_index * A_row_num; row_index < (col_index + 1) * A_row_num; row_index++) {
                double v = fabs(A[row_index]);
                // 1-norm of each column in A
                colnorms[col_index] += v;

                // compute cumsum of A
                cumsums[row_index + 1] = cumsums[row_index] + v;
            }
        }
    }

    /**
     * Compute cumulative sums and column norms of dense vector
     * @param A
     * @param A_col_num
     * @param colnorms
     * @param cumsums
     */
    inline void cumulative_sums(const double *A, int A_col_num, double *colnorms, double *cumsums) {
        cumsums[0] = 0.0;
        for (int col_index = 0; col_index < A_col_num; col_index++) {
            double v = fabs(A[col_index]);
            colnorms[col_index] = v;
            cumsums[col_index + 1] = cumsums[col_index] + v;
        }
    }

    /**
     * finds largest index of sorted (ascending) array a whose value is less than s
     * @param a
     * @param lb
     * @param ub
     * @param s
     * @return
     */
    inline int bin_search(double *a, int lb, int ub, double s) {
        int m;

        if (s < a[lb] || s > a[ub]) {
            cerr << "bin search: s out of bounds" << endl;
            exit(1);
        }

        while (lb < ub-1) {
            m = (lb + ub )/2;
            if (s < a[m]) {
                ub = m;
            } else {
                lb = m;
            }
        }

        return lb;
    }

    /**
     * https://stackoverflow.com/a/4609795/6620623
     * @tparam T
     * @param val
     * @return
     */
    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

public:

    WedgeNNSampling(const int QNum, const int PNum, const int d, double *QData, double *PData):QNum(QNum), PNum(PNum), d(d), QData(QData), PData(PData) {}

    ~WedgeNNSampling(){}

    void topk(const int k, const int s, const bool verify, const bool optimize, VectorElement *wedgeNNSampleResults) {

        double total_time = 0;

        Monitor timer;
        timer.start();

        double *col_norms_Q = new double[d]();
        double *col_norms_P = new double[d]();
        double *cumsums_Q = new double[d + 1]();
        double *cumsums_P = new double[PNum * d + 1]();
        int *lefts = new int[s]();
        int *rights = new int[s]();
        double *W = new double[d + 1]();
        int *bins = nullptr;
        int *centers = nullptr;

        if(optimize){
            bins = new int[d]();
        } else {
            centers = new int[s]();
        }

        VectorElement* verify_heap = nullptr;

        if(verify) {
            verify_heap = new VectorElement[s];
        }

        cumulative_sums(PData, PNum, d, col_norms_P, cumsums_P);

        for (int queryIndex = 0; queryIndex < QNum; queryIndex++) {

            double *qPtr = QData + queryIndex * d;

            cumulative_sums(qPtr, d, col_norms_Q, cumsums_Q);

            W[0] = 0.0;
            for (int l = 0; l < d; l++) {
                W[l + 1] = W[l] + col_norms_Q[l] * col_norms_P[l];
            }
            double tot_wt = W[d];

            if(optimize){
                std::fill(bins, bins + d, 0);
            }

            // choose centers
            for(int i = 0; i < s; i++) {
                double x = ((double)rand() / (double) RAND_MAX) * tot_wt;
                int e = bin_search(W, 0, d, x);
                if(optimize){
                    bins[e]++;
                } else {
                    centers[i] = e;
                }
            }

            // choose left
            if(optimize){
                int t = 0;
                for(int kk = 0; kk < d; kk++) {
                    for (int l = 0; l < bins[kk]; l++) {
                        double x = ((double) rand() / (double) RAND_MAX) * col_norms_Q[kk] + cumsums_Q[kk];
                        lefts[t] = bin_search(cumsums_Q, kk, kk + 1, x);
                        t++;
                    }
                }
            } else{
                for (int i = 0; i < s; i++) {
                    int kk = centers[i];
                    double x = ((double) rand() / (double) RAND_MAX) * col_norms_Q[kk] + cumsums_Q[kk];
                    lefts[i] = bin_search(cumsums_Q, kk, kk + 1, x);
                }
            }

            // choose right
            if(optimize){
                int t = 0;
                for(int jj = 0; jj < d; jj++) {
                    for (int l = 0; l < bins[jj]; l++) {
                        double x = ((double) rand() / (double) RAND_MAX) * col_norms_P[jj] + cumsums_P[jj * PNum];
                        rights[t] = bin_search(cumsums_P, jj * PNum, (jj + 1) * PNum, x);
                        t++;
                    }
                }
            } else {
                for (int i = 0; i < s; i++) {
                    int jj = centers[i];
                    double x = ((double) rand() / (double) RAND_MAX) * col_norms_P[jj] + cumsums_P[jj * PNum];
                    rights[i] = bin_search(cumsums_P, jj * PNum, (jj + 1) * PNum, x);
                }
            }

            // set output entry
            unordered_map<int, double> scores;

            for (int i = 0; i < s; i++) {
                int ii = lefts[i];
                int ll = rights[i] % PNum;
                scores[ll] += sgn(qPtr[lefts[i]]) * sgn(PData[rights[i]]);
            }

            VectorElement *heap = &wedgeNNSampleResults[queryIndex * k];

            if(verify) {
                EvalUtil::verified_retrieve(k, s, d, queryIndex, QData, PData, verify_heap, heap, scores);
            } else {
                EvalUtil::retrieve(k, heap, scores);
            }
        }

        timer.stop();

        total_time += timer.getElapsedTime();
        cout << "total time: " << total_time << " secs" << endl;

        delete[] col_norms_Q;
        delete[] col_norms_P;
        delete[] cumsums_Q;
        delete[] cumsums_P;
        delete[] lefts;
        delete[] rights;
        delete[] centers;
        delete[] bins;
        delete[] verify_heap;

    }
};
#endif //WEDGENNSAMPLING_H
