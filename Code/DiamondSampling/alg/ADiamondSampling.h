#ifndef DIAMOND_MIPS_T_H
#define DIAMOND_MIPS_T_H

#include "../util/Base.h"
#include "../util/CompareUtil.h"
#include "../util/Monitor.h"
#include "../util/RandomUtil.h"
#include "../structs/htree_t.h"
#include "../structs/VectorElement.h"
#include "../structs/htree_t.h"
#include "../structs/FastHeap.h"

/**
 *  Another DiamondSampling
 */
class ADiamondSampling {

private:

    int QNum, PNum, d;
    double *QData = nullptr;
    double *PData = nullptr;

    typedef htree_t<double> pool_t;
    vector<pool_t> P_cols;
    vector<int> g_cnt;
    vector<double> g_val;

    random_number_generator<> rng;

    void search(double *q, int budget, vector<int> &candidate) { // {{{

        htree_t<double> Ft(d), Forig(d);
        for(int t = 0; t < d; t++) {
            auto &P_col = P_cols[t];
            Ft[t] = fabs(q[t])*P_col.total_sum();
            Forig[t] = fabs(q[t]);
        }
        Ft.init_dense();
        for(int b = 0; b < budget; b++) {
            auto urnd = rng.uniform()*Ft.total_sum();
            int k = Ft.log_sample(urnd);
            auto &P_col = P_cols[k];
            urnd = rng.uniform()*P_col.total_sum();
            auto pIndex = P_col.log_sample(urnd);
            if(pIndex >= PNum) printf("sample error");
            urnd = rng.uniform()*Forig.total_sum();
            auto kPrime = Forig.log_sample(urnd);
            if(g_cnt[pIndex] == 0)
                candidate.push_back(pIndex);

            double *rowPtr1 = PData + pIndex * d;
            double *rowPtr2 = PData + kPrime * d;
            g_val[pIndex] += ((q[k]*rowPtr1[k]*q[kPrime]) > 0? 1 : -1) * rowPtr2[kPrime];
            g_cnt[pIndex] += 1;
        }

    }

    void search_true(double *q, const int k, const int budget, VectorElement *heap) { // {{{

        vector<int> candidate;
        search(q, budget, candidate);

        int heap_count = 0;
        for (auto &i : candidate) {
            double *pRowPtr = PData + i * d;
            double ip = std::inner_product(q, q + d, pRowPtr, 0.0);

            if (heap_count < k) {
                heap_enqueue(ip, i, heap, &heap_count);
            } else {
                if (ip > heap[0].data) {
                    heap_dequeue(heap, &heap_count);
                    heap_enqueue(ip, i, heap, &heap_count);
                }
            }
        }
//        std::sort(candidate.begin(), candidate.end(), CompareUtil::comparator<double>(g_val.data()));
        for (auto &i : candidate) {
            g_cnt[i] = 0;
//            g_val[i] = 0;
        }


    }

public:

    ADiamondSampling(const int QNum, const int PNum, const int d, double *QData, double *PData):QNum(QNum), PNum(PNum), d(d), QData(QData), PData(PData) {

        P_cols.resize(d, pool_t(PNum));
        g_cnt.resize(PNum);
        g_val.resize(PNum);

        for (int i = 0; i < d; i++) {
            auto &pool = P_cols[i];

            for (int j = 0; j < PNum; j++) {
                double *rowPtr = PData + j * d;
                pool[j] = fabs(rowPtr[i]);
            }
            pool.init_dense();
        }
    }

    ~ADiamondSampling(){}

    void topk(const int k, const int budget, VectorElement *sampledResults) {

        for (int queryIndex = 0; queryIndex < QNum; queryIndex++) {
            double *qPtr = QData + queryIndex * d;
            search_true(qPtr, k, budget, sampledResults + queryIndex * k);
        }
    }

};

#endif //DIAMOND_MIPS_T_H
