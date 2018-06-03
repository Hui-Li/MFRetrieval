#ifndef BUDGETEDMIP_H
#define BUDGETEDMIP_H

#include "../util/Base.h"
#include "../structs/FTree.h"
#include "../structs/arg_max.h"
#include "../util/Comparator.h"

class BudgetedMIP {

private:
    double *Q, *P;
    int QNum, PNum, d;
    typedef vector<arg_max> pool_t;
    vector<pool_t> pools;
    vector<size_t> g_cnt;
    vector<double> g_val;

    struct pool_view_t {
        const pool_t *p;
        bool increasing;
        size_t idx;

        pool_view_t() : p(NULL), increasing(true), idx(0) {}

        pool_view_t(const pool_t &p, bool increasing = true) : p(&p), increasing(increasing), idx(0) {}

        const arg_max &head() const { return (*this)[idx]; }

        bool empty() { return idx == p->size(); }

        void pop() { idx++; }

        const arg_max &operator[](size_t idx) const {
            if (increasing) return (*p)[idx];
            else return (*p)[p->size() - idx - 1];
        }
    };

    double search(double *q, size_t budget, vector<size_t> &candidate) { // {{{

        vector<pool_view_t> pool_views(d);
        Ftree_t<arg_max> Ft(d);
        for (size_t t = 0; t < d; t++) {
            pool_views[t] = pool_view_t(pools[t], q[t] < 0);
            Ft[t] = arg_max(t, pool_views[t].head().value * q[t]);
        }
        Ft.init_dense();
        budget = std::min(budget, pools.size() * pools[0].size() - 1);
        for (size_t b = 0; b < budget; b++) {
            auto t = Ft.val[1].idx;     // pool id
            auto v = Ft.val[1].value;   // entry value
            auto &cur_pool = pool_views[t];
            auto i = cur_pool.head().idx; // candidate id
            if (g_cnt[i] == 0)
                candidate.push_back(i);
            g_cnt[i] += 1;
            g_val[i] += v;

            cur_pool.pop();
            if (cur_pool.empty())
                Ft.set_value(t, arg_max());
            else
                Ft.set_value(t, arg_max(t, cur_pool.head().value * q[t]));
        }

    }

    void search_true(double *q, size_t budget, vector<size_t> &candidate) { // {{{

        candidate.clear();
        search(q, budget, candidate);
        for (auto &i : candidate) {
            double *pRowPtr = P + i * d;
            g_val[i] = std::inner_product(q, q + d, pRowPtr, 0.0);
        }
        std::sort(candidate.begin(), candidate.end(), comparator<double>(g_val.data()));
        for (auto &i : candidate) {
            g_cnt[i] = 0;
            g_val[i] = 0;
        }
    }

public:

    BudgetedMIP(const int QNum, const int PNum, const int d, double *Q, double *P): QNum(QNum), PNum(PNum), d(d), Q(Q), P(P) {
        pools.resize(d, pool_t(PNum));
        g_cnt.resize(PNum);
        g_val.resize(PNum);

        for (size_t i = 0; i < d; i++) {
            auto &pool = pools[i];
            for (size_t j = 0; j < PNum; j++) {
                double *rowPtr = P + j * d;
                pool[j] = arg_max(j, rowPtr[i]);
            }
            std::sort(pool.begin(), pool.end()); // increasing order
        }
    }

    void topk(const int k, const int budget, vector<vector<size_t> > &candidates) {

        for(size_t qIndex = 0; qIndex < QNum; qIndex++) {
            double *QPtr = Q + qIndex * d;
            vector<size_t> &candidate = candidates[qIndex];
            search_true(QPtr, budget, candidate);
        }

    }
};

#endif //BUDGETEDMIP_H
