#ifndef DARYBASESCHEMA_H
#define DARYBASESCHEMA_H

#include "../util/Base.h"
#include "../structs/Matrix.h"
#include "../util/CompareUtil.h"
#include "../structs/VectorElement.h"
#include "../structs/FastHeap.h"

class DaryBaseSchema{

public:

    DaryBaseSchema(Matrix &Q, Matrix &P, const int d) : Q(Q), P(P), d(d) {
        original_dim = Q.colNum;
        D = 2 * d + 1;
        p = original_dim * D;
    }

    inline void data_prepare(){

#ifdef DEBUG
        cout << "d:" << d << endl;
        cout << "D:" << D << endl;
        cout << "p:" << p << endl;
#endif

        base.resize(D);

        // {1, 0.99, 0.98, ..., 0, -0.01, ...-1}

        for (int i = 0; i <= d; i++) {
            base[i] = (d - i) / (d + 0.0);
        }

        for (int i = d + 1; i <= 2 * d; i++) {
            base[i] = (d - i) / (d + 0.0);
        }

        // {-1, -0.99, -0.98, ..., 0, 0.01, ...1}

//        for (int i = 0; i <= d; i++) {
//            base[i] = -(d - i) / (d + 0.0);
//        }
//
//        for (int i = d + 1; i <= 2 * d; i++) {
//            base[i] = -(d - i) / (d + 0.0);
//        }

#ifdef DEBUG
        cout << "base" << endl;
        for(auto b:base){
            cout << b << " ";
        }
        cout << endl;
        cout << "-----------------------" << endl;
#endif

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
            tess_vector_dary(Q[query_id], q_a_row_ptr);

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
    vector<double> base;
    int p;
    int D, d;
    int original_dim;

    inline void permutate(Matrix &matrix, Matrix &sparse_matrix, vector<set<int> > &inverted_indices) {
        inverted_indices.resize(matrix.rowNum);
        double *a_row_ptr = new double[original_dim];
        for (int row_index = 0; row_index < matrix.rowNum; row_index++) {
            tess_vector_dary(matrix[row_index], a_row_ptr);
            permutate_vec(inverted_indices[row_index], matrix[row_index], a_row_ptr, sparse_matrix[row_index]);
        }
        delete[] a_row_ptr;
    }

    int count = 0;

    inline void permutate_vec(set<int> &inverted_index, double *original_row_ptr, double *a_row_ptr, double *sparse_row_ptr){


#ifdef DEBUG
        cout << "a" << endl;
        for (int i = 0; i < original_dim; i++) {
            cout << a_row_ptr[i] << " ";
        }
        cout << endl;

        cout << "original" << endl;
        for (int i = 0; i < original_dim; i++) {
            cout << original_row_ptr[i] << " ";
        }
        cout << endl;

        cout << "sparse" << endl;
#endif

        for (int i = 0; i < p; i++) {
            int original_index = i / D;
            int remainder = i % D;

            if(a_row_ptr[original_index] == base[remainder]){
                sparse_row_ptr[i] = original_row_ptr[original_index];
            } else {
                sparse_row_ptr[i] = 0;
            }

            if(sparse_row_ptr[i]!=0){
                inverted_index.insert(i);
            }

#ifdef DEBUG
            cout << remainder << ":" << a_row_ptr[original_index] << "==" << base[remainder] << "->" << sparse_row_ptr[i] << "|";
#endif

        }

#ifdef DEBUG
        cout << endl;

        for (int i = 0; i < p; i++) {
            cout << sparse_row_ptr[i] << " ";
        }

        cout << endl;

        cout << "-------------------" << endl;

        count ++;
        if(count % 10==0) {
            exit(1);
        }
#endif

    }

    inline void tess_vector_dary(double *row_ptr, double *a_ptr) {

        for (int i = 0; i < original_dim; i++) {
            double upper = std::ceil(d * row_ptr[i]);
            double lower = std::floor(d * row_ptr[i]);
            double a_plus = std::abs(d * row_ptr[i] - upper);
            double a_minus = std::abs(d * row_ptr[i] - lower);
            if (a_plus <= a_minus) {
                a_ptr[i] = upper / d;
            } else {
                a_ptr[i] = lower / d;
            }
        }

//        double norm = 0;
//        for (int i = 0; i < original_dim; i++) {
//            norm += a_ptr[i] * a_ptr[i];
//        }
//
//        norm = sqrt(norm);
//
//        for (int i = 0; i < original_dim; i++) {
//            a_ptr[i] /= norm;
//        }

    }

};

#endif //DARYBASESCHEMA_H
