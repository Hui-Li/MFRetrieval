#ifndef TERNARYBASESCHEMA_H
#define TERNARYBASESCHEMA_H

#include "../util/Base.h"
#include "../structs/Matrix.h"
#include "../util/CompareUtil.h"
#include "../structs/VectorElement.h"
#include "../structs/FastHeap.h"

class TernaryBaseSchema{

public:

    TernaryBaseSchema(Matrix &Q, Matrix &P, const int p) : Q(Q), P(P), p(p) {
        original_dim = Q.colNum;
    }

    inline void data_prepare(){

        sparse_P.init(P.rowNum, p);
        permutate(P, sparse_P, P_inverted_index);

#ifdef DEBUG
        Matrix sparse_Q;
        sparse_Q.init(Q.rowNum, p);
        permutate(Q, sparse_Q, Q_inverted_index);

        string q_inv = "q_inv.txt";
        string p_inv = "p_inv.txt";

        ofstream q_file(q_inv.c_str());
        for (int i = 0; i < Q_inverted_index.size(); i++) {
            set<int> &s = Q_inverted_index[i];
            for (const int &s_e:s) {
                q_file << s_e << " ";
            }
            q_file << endl;
        }
        q_file.close();

        ofstream p_file(p_inv.c_str());
        for (int i = 0; i < P_inverted_index.size(); i++) {
            set<int> &s = P_inverted_index[i];
            for (const int &s_e:s) {
                p_file << s_e << " ";
            }
            p_file << endl;
        }
        p_file.close();

#endif

    }

    inline void query_topk(const int k, VectorElement *results){

        vector<int> common_data;
        int heapCount = 0;

        double *q_a_row_ptr = new double[original_dim];
        double *sparse_q_ptr = new double[p];

        for (int query_id = 0; query_id < Q.rowNum; query_id++) {
            if(query_id % 1000 == 0){
                cout << query_id << " of " << Q.rowNum << endl;
            }

            set<int> q_inverted_index;
            tess_vector(Q[query_id], q_a_row_ptr);

            permutate_vec(q_inverted_index, Q[query_id], q_a_row_ptr, sparse_q_ptr);

            VectorElement *heap = &results[query_id * k];
            heapCount = 0;

            for (int p_row_id = 0; p_row_id < P.rowNum; p_row_id++) {
                common_data.clear();
                set<int> &p_set = P_inverted_index[p_row_id];

                set_intersection(q_inverted_index.begin(), q_inverted_index.end(), p_set.begin(), p_set.end(),
                                 std::back_inserter(common_data));

                if (common_data.size() > 0) {

                    double inner_product = 0;
                    double *p_row_ptr = sparse_P[p_row_id];
                    for(auto dim:common_data){
                        inner_product += sparse_q_ptr[dim] * p_row_ptr[dim];
                    }

                    if (heapCount < k) {
                        heap_enqueue(inner_product, p_row_id, heap, &heapCount);
                    } else if (inner_product > heap[0].data) {
                        heap_dequeue(heap, &heapCount);
                        heap_enqueue(inner_product, p_row_id, heap, &heapCount);
                    }
                }
            }
        }

        delete[] q_a_row_ptr;
        delete[] sparse_q_ptr;
    }

private:

    Matrix &Q, &P;
    Matrix sparse_P;
    vector<set<int> > Q_inverted_index;
    vector<set<int> > P_inverted_index;
    int p;
    int original_dim;

    inline void permutate(Matrix &matrix, Matrix &sparse_matrix, vector<set<int> > &inverted_indices) {
        inverted_indices.resize(matrix.rowNum);
        double *a_row_ptr = new double[original_dim];
        for (int row_index = 0; row_index < matrix.rowNum; row_index++) {
            tess_vector(matrix[row_index], a_row_ptr);
            permutate_vec(inverted_indices[row_index], matrix[row_index], a_row_ptr, sparse_matrix[row_index]);
        }
        delete[] a_row_ptr;
    }

    //    int count = 0;

    inline void permutate_vec(set<int> &inverted_index, double *original_row_ptr, double *a_row_ptr, double *sparse_row_ptr){

//        for (int i = 0; i < original_dim; i++) {
//            cout << original_row_ptr[i] << " ";
//        }
//        cout << endl;

        for (int i = 0; i < p; i++) {
            int original_index = i / 3;

            if (((i % 3) == 0) && (a_row_ptr[original_index] > 0)) { // a_row_ptr[original_index] == 1
                sparse_row_ptr[i] = original_row_ptr[original_index];
            } else if (((i % 3) == 1) && (a_row_ptr[original_index] == 0)) {
                sparse_row_ptr[i] = original_row_ptr[original_index];
            } else if (((i % 3) == 2) && (a_row_ptr[original_index] < 0)) { // a_row_ptr[original_index] == -1
                sparse_row_ptr[i] = original_row_ptr[original_index];
            } else {
                sparse_row_ptr[i] = 0;
            }

            if(sparse_row_ptr[i]!=0){
                inverted_index.insert(i);
            }

//            cout << sparse_row_ptr[i] << " ";
        }
//        cout << endl;
//
//        cout << "-------------------" << endl;
//
//        count ++;
//        if(count % 10==0) {
//            exit(1);
//        }
    }

    inline void tess_vector(double *original_row_ptr, double *a_ptr) {

        vector<pair<int, double> > row_vec(original_dim);

        for (int i = 0; i < original_dim; i++) {
            row_vec[i].first = i;
            row_vec[i].second = std::abs(*(original_row_ptr + i));
        }

        sort(row_vec.begin(), row_vec.end(), CompareUtil::pairGreaterCompare);

        double cum_sum = 0;
        double max_value = -1;
        int max_index = -1;

        set<int> index_set;

        for (int i = 0; i < original_dim; i++) {
            cum_sum += row_vec[i].second;
            double curr_value = cum_sum/sqrt(i+1);
            if (curr_value > max_value) {
                max_index = i;
                max_value = curr_value;
            }

            index_set.insert(max_index);
        }

        double s = sqrt(index_set.size());

#ifdef DEBUG
        double norm = 0;
#endif
        for (int i = 0; i < original_dim; i++) {
            if (index_set.find(i) != index_set.end()) {
                a_ptr[i] = sign(original_row_ptr[i]) / s;
            } else {
                a_ptr[i] = 0;
            }
#ifdef DEBUG
            norm += a_ptr[i] * a_ptr[i];
#endif
        }

#ifdef DEBUG
        // norm should be one
        assert(sqrt(norm)==1);
#endif
    }

};

#endif //TERNARYBASESCHEMA_H
