#ifndef EVALUTIL_H
#define EVALUTIL_H


// average recall, see paper:
// EFANNA : An Extremely Fast Approximate Nearest Neighbor Search Algorithm Based on kNN Graph
// Section 4.2.1
void avg_recall(const int QNum, const int k, string groundTruthFilePath, VectorElement *diamondSampleResults) {

    vector<pair<unordered_set<int>, unordered_set<double> > > groundTruth(QNum, std::make_pair(unordered_set<int>(), unordered_set<double>()));

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
                pair<unordered_set<int>, unordered_set<double> > &truth = groundTruth[qIndex];
                VectorElement *result = diamondSampleResults + qIndex * k;

                for (int i = 0; i < k; i++) {
                    if ((truth.first.find(result[i].id) != truth.first.end()) ||
                        (truth.second.find(result[i].data) != truth.second.end())) {
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

    cout << "avg hit: " << total_count / (QNum * 1.0) << endl;
    double recall = total_count / (k * QNum * 1.0);

    cout << "avg recall: " << recall << endl;
}


#endif //EVALUTIL_H
