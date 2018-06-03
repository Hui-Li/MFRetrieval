#ifndef EVALUTIL_H
#define EVALUTIL_H

#include "Base.h"
#include "FileUtil.h"
#include "../structs/VectorElement.h"

namespace EvalUtil {

    // average recall, see paper:
    // EFANNA : An Extremely Fast Approximate Nearest Neighbor Search Algorithm Based on kNN Graph
    void avg_recall(const int QNum, const int k, string groundTruthFilePath, vector<vector<size_t> > &results) {
        vector<unordered_set<int> > groundTruth(QNum);
        FileUtil::readGroundtruth(groundTruth, groundTruthFilePath);

        int threadNum = std::thread::hardware_concurrency();

        vector<int> count(threadNum, 0);

        vector<std::thread> exec_threads(threadNum);
        int workload = QNum / threadNum;

        for (int thread_index = 0; thread_index < threadNum; thread_index++) {

            exec_threads[thread_index] = std::thread(std::bind([&](const int thread_index) {
                int start = thread_index * workload;
                int end = (thread_index + 1) * workload;
                end = (end > QNum) ? QNum : end;

                for (int qIndex = start; qIndex < end; qIndex++) {
                    unordered_set<int> &truth = groundTruth[qIndex];
                    vector<size_t> &result = results[qIndex];

                    for (int i = 0; i < k; i++) {
                        if (truth.find(result[i]) != truth.end()) {
                            count[thread_index]++;
                        }
                    }
                }

            }, thread_index));
        }

        for (int thread_index = 0; thread_index < threadNum; thread_index++) {
            exec_threads[thread_index].join();
        }

        int total_count = std::accumulate(count.begin(), count.end(), 0);

        double recall = total_count / (k * QNum * 1.0);

        cout << "avg recall: " << recall << endl;
    }

    inline void verified_retrieve(const int k, const int search_k, const int d, const int qIndex, double *QData, double *PData,
                      VectorElement *large_heap, VectorElement *small_heap,
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

    inline void retrieve(const int k, VectorElement *heap,
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
}

#endif //EVALUTIL_H
