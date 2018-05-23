#ifndef DIAMONDSAMPLING_H
#define DIAMONDSAMPLING_H

#include <random>
#include "../structs/VectorElement.h"
#include "../util/Base.h"
#include "../util/FileUtil.h"
#include "../util/AliasSamplingGenerator.h"
#include "../structs/FastHeap.h"
#include "../util/Monitor.h"

#define DEBUG

class DiamondSampling{

private:

    int QNum, PNum, d;
    double *QData = nullptr;
    double *PData = nullptr;

    inline void verified_retrieve(const int k, const int search_k, const int qIndex, VectorElement *large_heap, VectorElement *small_heap,
                         unordered_map<int, double> &sampled_score) {

        int large_heap_count = 0;

        for (auto ptr = sampled_score.begin(); ptr != sampled_score.end(); ptr++) {
            if (large_heap_count < search_k) {
                heap_enqueue(ptr->second, ptr->first, large_heap, &large_heap_count);
            } else {
                if (ptr->second > large_heap[0].data) {
                    heap_dequeue(large_heap, &large_heap_count);
                    heap_enqueue(ptr->second, ptr->first, large_heap, &large_heap_count);
                }
            }
        }

        int small_heap_count = 0;
        for (int i = 0; i < large_heap_count; i++) {

            double *qRowPtr = QData + qIndex * d;
            int pIndex = large_heap[i].id;
            double *pRowPtr = PData + pIndex * d;
            double ip = std::inner_product(qRowPtr, qRowPtr + d, pRowPtr, 0.0);

            if (small_heap_count < k) {
                heap_enqueue(ip, pIndex, small_heap, &small_heap_count);
            } else {
                if (ip > small_heap[0].data) {
                    heap_dequeue(small_heap, &small_heap_count);
                    heap_enqueue(ip, pIndex, small_heap, &small_heap_count);
                }
            }

        }

    }

    inline void retrieve(const int k, const int qIndex, VectorElement *heap,
                                  unordered_map<int, double> &sampled_score) {
        int heap_count = 0;

        for (auto ptr = sampled_score.begin(); ptr != sampled_score.end(); ptr++) {
            if (heap_count < k) {
                heap_enqueue(ptr->second, ptr->first, heap, &heap_count);
            } else {
                if (ptr->second > heap[0].data) {
                    heap_dequeue(heap, &heap_count);
                    heap_enqueue(ptr->second, ptr->first, heap, &heap_count);
                }
            }
        }
    }

public:

    DiamondSampling(const int QNum, const int PNum, const int d, double *QData, double *PData):QNum(QNum), PNum(PNum), d(d), QData(QData), PData(PData) {
    }

    ~DiamondSampling(){}

    void topk(const int k, const int search_k, const int s, const bool verify, VectorElement *sampledResults) {

        double total_time = 0;

        Monitor timer;

        timer.start();

        vector<double> edges_weights(d);
        vector<double> q_probability(d);
        vector<vector<double> > p_probabilities(d, vector<double>(PNum));
        vector<double> P_column_sum(d, 0);
        unordered_map<int, double> sampled_score;
        VectorElement* large_heap = nullptr;
        if(verify) {
            large_heap = new VectorElement[search_k];
        }

        for (int pIndex = 0; pIndex < PNum; pIndex++) {
            double *rowPtr = PData + pIndex * d;
            for (int col = 0; col < d; col++) {
                double value = std::abs(rowPtr[col]);
                p_probabilities[col][pIndex] = value;
                P_column_sum[col] += value;
            }
        }

        vector<AliasSamplingGenerator *> samplers2(d, nullptr);
        for (int col = 0; col < d; col++) {
            vector<double> &p_probability = p_probabilities[col];
            samplers2[col] = new AliasSamplingGenerator(p_probability.size(), p_probability.data());
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1.0);

        timer.stop();
        double prepare_time = timer.getElapsedTime();
        cout << "prepare_time: " << prepare_time << " secs" << endl;

        double prepare_time2 = 0;
        double sample_time = 0;
        double retrieve_time = 0;

        for (int queryIndex = 0; queryIndex < QNum; queryIndex++) {

//            if(queryIndex % 1000 == 0) {
//                cout << queryIndex << " of " << QNum << endl;
//            }

            timer.start();

            sampled_score.clear();

            double *qPtr = QData + queryIndex * d;
            double q_Norm = 0;
            for (int col = 0; col < d; col++) {
                q_probability[col] = std::abs(qPtr[col]);
                q_Norm += q_probability[col];
            }

            for (int col = 0; col < d; col++) {
                edges_weights[col] = q_probability[col] * q_Norm * P_column_sum[col];
            }

            timer.stop();
            prepare_time2 += timer.getElapsedTime();

            timer.start();
            AliasSamplingGenerator sampler1(edges_weights.size(), edges_weights.data());
            AliasSamplingGenerator sampler3(q_probability.size(), q_probability.data());

            for (int i = 0; i < s; i++) {
                // line 6 in Alg 3
                int k = sampler1.sample(dis(gen), dis(gen));

                // line 7 in Alg 3
                AliasSamplingGenerator &sampler2 = * (samplers2[k]);
                int pIndex = sampler2.sample(dis(gen), dis(gen));

                // line 8 in Alg3
                int kPrime = sampler3.sample(dis(gen), dis(gen));

                double *pRowPtr = PData + pIndex * d;

                double sgn = qPtr[k] * pRowPtr[k] * qPtr[kPrime];
                sgn = (sgn >= 0) ? 1 : -1;
                if(sampled_score.find(pIndex) == sampled_score.end()){
                    sampled_score[pIndex] = sgn * pRowPtr[kPrime];
                } else {
                    sampled_score[pIndex] += sgn * pRowPtr[kPrime];
                }
            }

            timer.stop();
            sample_time += timer.getElapsedTime();

            timer.start();

            if(verify) {
                verified_retrieve(k, search_k, queryIndex, large_heap, sampledResults + k * queryIndex, sampled_score);
            } else {
                retrieve(k, queryIndex, sampledResults + k * queryIndex, sampled_score);
            }
            timer.stop();
            retrieve_time += timer.getElapsedTime();
        }

        if(verify) {
            delete[] large_heap;
        }

        for (int i = 0; i < samplers2.size(); i++) {
            delete samplers2[i];
        }
        cout << "prepare_time2: " << prepare_time2 << " secs" << endl;
        cout << "sample_time: " << sample_time << " secs" << endl;
        cout << "retrieve_time: " << retrieve_time << " secs" << endl;

        cout << "total time: " << prepare_time + prepare_time2 + sample_time + retrieve_time << " secs" << endl;

    }
};
#endif //DIAMONDSAMPLING_H
