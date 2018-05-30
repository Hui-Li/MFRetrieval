#ifndef WEDGESAMPLING_H
#define WEDGESAMPLING_H

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


class WedgeSampling{

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

    inline void compute_prob(Node *node){
        for (int col = 0; col < node->next_level_size; col++) {
            Edge *edge = node->to_next_level[col];
            node->transistion[col] = edge->weight * edge->to->W / edge->from->W;
        }

        node->clear_generator();
    }

    inline void compute_probs(const int level_size, Node *level){
        for (int i = 0; i < level_size; i++) {
            compute_prob(level + i);
        }
    }

    inline int sample(Node *start_node){

        // first -> second
        int target_edge_index = DiscreteRandomGenerator::DiscreteRandom(start_node);
        Node *node_second_level = start_node->to_next_level[target_edge_index]->to;

        // second -> third
        target_edge_index = DiscreteRandomGenerator::DiscreteRandom(node_second_level);

        return target_edge_index;
    }

public:

     WedgeSampling(const int QNum, const int PNum, const int d, double *QData, double *PData):QNum(QNum), PNum(PNum), d(d), QData(QData), PData(PData) {
    }

    ~WedgeSampling(){
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
        if(verify){
            verify_heap = new VectorElement[s]();
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

        // Compute W and Probability for second level -> third level
        compute_Ws(second_level_size, P_second_level);
        compute_Ws(second_level_size, N_second_level);

        compute_probs(second_level_size, P_second_level);
        compute_probs(second_level_size, N_second_level);

        timer.stop();
        cout << "prepare time: " << timer.getElapsedTime() << " secs" << endl;
        total_time += timer.getElapsedTime();

        timer.start();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1.0);

        vector<double> weights(2);

#ifdef DEBUG
        double set_edge_time = 0;
        double compute_w_time = 0;
        double clear_map_time = 0;
        double nonuniformed_sample_time = 0;
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
            compute_probs(first_level_size, P_first_level);
            compute_probs(first_level_size, N_first_level);

            double positive_weight = P_first_level[0].W;
            double negative_weight = N_first_level[1].W;

#ifdef DEBUG
            timer2.stop();
            compute_w_time += timer2.getElapsedTime();
#endif

#ifdef DEBUG
            timer2.start();
#endif
            unordered_map<int, double> scores;
//            scores.clear();
#ifdef DEBUG
            timer2.stop();
            clear_map_time += timer2.getElapsedTime();
#endif

            for (int i = 0; i < s; i++) {
#ifdef DEBUG
                timer2.start();
#endif

                int pos_sample = sample(P_first_level);
                int neg_sample = sample(N_first_level);

#ifdef DEBUG
                timer2.stop();
                nonuniformed_sample_time += timer2.getElapsedTime();

                timer2.start();
#endif
                scores[pos_sample] += positive_weight;
                scores[neg_sample] -= negative_weight;

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
        cout << "clear_map_time: " << clear_map_time << " secs" << endl;
        cout << "nonuniformed_sample_time: " << nonuniformed_sample_time << " secs" << endl;
        cout << "search_time: " << search_time << " secs" << endl;
        cout << "heap_time: " << heap_time << " secs" << endl;
        cout << "------------------------" << endl;
#endif
    }
};
#endif //WEDGESAMPLING_H
