#ifndef WEDGEPRESAMPLING_H
#define WEDGEPRESAMPLING_H

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
#include "../util/EvalUtil.h"

class WedgePreSampling{

private:

    int QNum, PNum, d, s, first_level_size, second_level_size, third_level_size;
    double *QData = nullptr;
    double *PData = nullptr;
    double *P_T = nullptr;
    Node *P_first_level = nullptr;
    Node *N_first_level = nullptr;
    Node *P_second_level = nullptr;
    Node *N_second_level = nullptr;
    Node *third_level = nullptr;

    vector<vector<int> > positive_pre_sampled_results, negative_pre_sampled_results;

    vector<double> p_a, n_a;

    inline void set_edge(Node *from_node, Node *to_node, const int edge_index, const double weight){
        if(from_node->next_level_size > 0) {
            Edge *from_to_edge = from_node->to_next_level[edge_index];
            from_to_edge->from = from_node;
            from_to_edge->to = to_node;
            from_to_edge->weight = weight;
        }
    }

    inline void compute_W(Node *node){
        double W = 0;
        for (int col = 0; col < node->next_level_size; col++) {
            Edge *edge = node->to_next_level[col];
            W+= edge->weight * edge->to->W;
        }
        node->W = W;
    }

    inline void compute_Ws(const int level_size, Node *level){
        for (int i = 0; i < level_size; i++) {
            compute_W(level + i);
        }
    }

    inline void pre_sample(){

        p_a.resize(d, 0);
        n_a.resize(d, 0);

        vector<double> p_weight(PNum, 0), n_weight(PNum, 0);
        positive_pre_sampled_results.resize(d, vector<int>(s));
        negative_pre_sampled_results.resize(d, vector<int>(s));

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1.0);

        for (int j = 0; j < d; j++) {

            p_weight.resize(PNum, 0);
            n_weight.resize(PNum, 0);

            double *rowPtr = P_T + j * PNum;

            for (int pIndex = 0; pIndex < PNum; pIndex++) {
                if (rowPtr[pIndex] > 0) {
                    p_a[j] += rowPtr[pIndex];
                    p_weight[pIndex] = rowPtr[j];
                } else {
                    n_a[j] -= rowPtr[pIndex];
                    n_weight[pIndex] = -rowPtr[j];
                }
            }

            vector<int> &p_pre_samples = positive_pre_sampled_results[j];
            vector<int> &n_pre_samples = negative_pre_sampled_results[j];

            AliasSamplingGenerator p_generator(p_weight.size(), p_weight.data());
            AliasSamplingGenerator n_generator(p_weight.size(), p_weight.data());

            for (int i = 0; i < s; i++) {
                p_pre_samples[i] = p_generator.sample(dis(gen), dis(gen));
                n_pre_samples[i] = n_generator.sample(dis(gen), dis(gen));
            }
        }

    }

public:

     WedgePreSampling(const int QNum, const int PNum, const int d, double *QData, double *PData):QNum(QNum), PNum(PNum), d(d), QData(QData), PData(PData) {
    }

    ~WedgePreSampling(){
        delete[] P_T;
        delete[] P_first_level;
        delete[] N_first_level;
        delete[] P_second_level;
        delete[] N_second_level;
        delete[] third_level;
    }

    void topk(const int k, const int s, const bool verify, VectorElement *wedgeSampleResults) {
        this->s = s;
        double total_time = 0;

        Monitor timer;

        timer.start();

        first_level_size = 1;
        second_level_size = d;
        third_level_size = PNum;

        VectorElement* verify_heap = nullptr;

        if(verify) {
            verify_heap = new VectorElement[s];
        }

        P_T = new double[d * PNum]; // d * PNum
        Calculator::transpose(PData, P_T, PNum, d);

        P_first_level = new Node[first_level_size];
        N_first_level = new Node[first_level_size];
        P_second_level = new Node[second_level_size];
        N_second_level = new Node[second_level_size];
        third_level = new Node[third_level_size];

        for (int i = 0; i < first_level_size; i++) {
            (P_first_level + i)->init(second_level_size);
            (N_first_level + i)->init(second_level_size);

            for (int j = 0; j < second_level_size; j++) {
                P_first_level[i].to_next_level[j] = new Edge();
                N_first_level[i].to_next_level[j] = new Edge();
            }
        }

        for (int i = 0; i < second_level_size; i++) {
            (P_second_level + i)->init(third_level_size);
            (N_second_level + i)->init(third_level_size);

            for (int j = 0; j < third_level_size; j++) {
                P_second_level[i].to_next_level[j] = new Edge();
                N_second_level[i].to_next_level[j] = new Edge();
            }
        }

        // set edges and edge weight w(e) from second level to third level
        for(int second_level_index = 0; second_level_index < second_level_size; second_level_index ++) {

            for (int third_level_index = 0; third_level_index < third_level_size; third_level_index++) {
                double value = P_T[second_level_index * PNum + third_level_index];
                if(value > 0) { // +node
                    set_edge(P_second_level + second_level_index, third_level + third_level_index, third_level_index, value);
                    set_edge(N_second_level + second_level_index, third_level + third_level_index, third_level_index, 0);
                } else {        // -node
                    set_edge(N_second_level + second_level_index, third_level + third_level_index, third_level_index, -value);
                    set_edge(P_second_level + second_level_index, third_level + third_level_index, third_level_index, 0);
                }
            }
        }

        // Compute W values for second level
        compute_Ws(second_level_size, P_second_level);
        compute_Ws(second_level_size, N_second_level);

        timer.stop();
        cout << "prepare time: " << timer.getElapsedTime() << " secs" << endl;
        total_time += timer.getElapsedTime();

        timer.start();
        pre_sample();
        timer.stop();
        cout << "presample time: " << timer.getElapsedTime() << " secs" << endl;
        total_time += timer.getElapsedTime();

        timer.start();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1.0);

        vector<double> weights(2);

#ifdef DEBUG
        double set_edge_time = 0;
        double compute_w_time = 0;
        double sum_time = 0;
        double nonuniformed_sample_time = 0;
        double round_time = 0;
        double search_time = 0;
        double heap_time = 0;
        Monitor timer2;
#endif

        for (int queryIndex = 0; queryIndex < QNum; queryIndex++) {

            double *qPtr = QData + queryIndex * d;

#ifdef DEBUG
            timer2.start();
#endif

            // set edges and edge weight w(e) from first level to second level
            for (int second_level_index = 0; second_level_index < second_level_size; second_level_index++) {
                double value = qPtr[second_level_index];
                if (value > 0) { // +node
                    // ++
                    set_edge(P_first_level, P_second_level + second_level_index, second_level_index,
                             value);
                    // --
                    set_edge(N_first_level, N_second_level + second_level_index, second_level_index,
                             value);
                } else {        // -node
                    // +-
                    set_edge(P_first_level, N_second_level + second_level_index, second_level_index,
                             -value);
                    // -+
                    set_edge(N_first_level, P_second_level + second_level_index, second_level_index,
                             -value);
                }
            }

#ifdef DEBUG
            timer2.stop();

            set_edge_time += timer2.getElapsedTime();

            timer2.start();
#endif

            // Compute W values for the first level (used together with sampled results for ranking)
            compute_Ws(first_level_size, P_first_level);
            compute_Ws(first_level_size, N_first_level);

#ifdef DEBUG
            timer2.stop();
            compute_w_time += timer2.getElapsedTime();

            timer2.start();
#endif

            double positive_weight = P_first_level[0].W;
            double negative_weight = N_first_level[1].W;

            vector<double> positive_preprocessed_query(d, 0);
            vector<double> negative_preprocessed_query(d, 0);

            unordered_map<int, double> scores;

            for (int i = 0; i < d; i++) {
                if (qPtr[i] > 0) {
                    positive_preprocessed_query[i] = qPtr[i] * p_a[i];
//                    positive_sum += positive_preprocessed_query[i];
                } else {
                    negative_preprocessed_query[i] = -qPtr[i] * n_a[i];
//                    negative_sum += negative_preprocessed_query[i];
                }
            }

            // this is faster!
            double positive_sum = std::accumulate(positive_preprocessed_query.begin(), positive_preprocessed_query.end(), 0.0);;
            double negative_sum = std::accumulate(negative_preprocessed_query.begin(), negative_preprocessed_query.end(), 0.0);;

#ifdef DEBUG
            timer2.stop();
            sum_time += timer2.getElapsedTime();
#endif

            for (int i = 0; i < d; i++) {

                if(positive_preprocessed_query[i]==0){
                    continue;
                }

#ifdef DEBUG
                timer2.start();
#endif

                double value = positive_preprocessed_query[i] / positive_sum * s;

                int floor_value = std::floor(value);
                int ceil_value = std::ceil(value);

                weights[0] = ceil_value - value;
                weights[1] = value - floor_value;

#ifdef DEBUG
                timer2.stop();
                round_time += timer2.getElapsedTime();

                timer2.start();
#endif

                AliasSamplingGenerator p_generator(weights.size(), weights.data());

                int positive_sample_size = p_generator.sample(dis(gen), dis(gen)) == 0 ? floor_value : ceil_value;

#ifdef DEBUG
                timer2.stop();
                nonuniformed_sample_time += timer2.getElapsedTime();

                timer2.start();
#endif

                vector<int> &p_pre_samples = positive_pre_sampled_results[i];

                for (int sample_index = 0; sample_index < positive_sample_size; sample_index++) {
                    int sampled_target = p_pre_samples[sample_index];
                    scores[sampled_target] += positive_weight;
                }

#ifdef DEBUG
                timer2.stop();
                search_time += timer2.getElapsedTime();
#endif

            }

            for (int i = 0; i < d; i++) {

                if(negative_preprocessed_query[i]==0){
                    continue;
                }

#ifdef DEBUG
                timer2.start();
#endif

                double value = s * negative_preprocessed_query[i] / negative_sum;

                int floor_value = std::floor(value);
                int ceil_value = std::ceil(value);

                weights[0] = ceil_value - value;
                weights[1] = value - floor_value;

#ifdef DEBUG
                timer2.stop();
                round_time += timer2.getElapsedTime();

                timer2.start();
#endif

                AliasSamplingGenerator n_generator(weights.size(), weights.data());

#ifdef DEBUG
                timer2.stop();
                nonuniformed_sample_time += timer2.getElapsedTime();
#endif

                int negative_sample_size = n_generator.sample(dis(gen), dis(gen)) == 0 ? floor_value : ceil_value;

#ifdef DEBUG
                timer2.start();
#endif

                vector<int> &n_pre_samples = negative_pre_sampled_results[i];

                for (int sample_index = 0; sample_index < negative_sample_size; sample_index++) {

                    int sampled_target = n_pre_samples[sample_index];
                    scores[sampled_target] -= negative_weight;

                }

#ifdef DEBUG
                timer2.stop();
                search_time += timer2.getElapsedTime();
#endif

            }

#ifdef DEBUG
            timer2.start();
#endif

            VectorElement *heap = &wedgeSampleResults[queryIndex * k];

            if(verify) {
                EvalUtil::verified_retrieve(k, s, d, queryIndex, QData, PData, verify_heap, heap, scores);
            } else {
                EvalUtil::retrieve(k, heap, scores);
            }

#ifdef DEBUG
            timer2.stop();
            heap_time += timer2.getElapsedTime();
#endif

        }

        timer.stop();
        cout << "retrieve time: " << timer.getElapsedTime() << " secs" << endl;
        total_time += timer.getElapsedTime();
        cout << "total time: " << total_time << " secs" << endl;

        delete[] verify_heap;
#ifdef DEBUG
        cout << "---- detailed cost ----" << endl;
        cout << "set_edge_time: " << set_edge_time << " secs" << endl;
        cout << "compute_w_time: " << compute_w_time << " secs" << endl;
        cout << "sum_time: " << sum_time << " secs" << endl;
        cout << "round_time: " << round_time << " secs" << endl;
        cout << "nonuniformed_sample_time: " << nonuniformed_sample_time << " secs" << endl;
        cout << "search_time: " << search_time << " secs" << endl;
        cout << "heap_time: " << heap_time << " secs" << endl;
        cout << "------------------------" << endl;
#endif
    }
};
#endif //WEDGEPRESAMPLING_H
