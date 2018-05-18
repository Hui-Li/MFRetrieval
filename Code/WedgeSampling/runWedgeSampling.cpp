#include "util/Base.h"
#include "alg/WedgeSampling.h"

/**
 * Naive solution.
 */
void naive(const int threadNum, const int qNum, const int pNum, const int d, const int k, const double *q, const double *p,
           VectorElement *results) {
    Monitor timer;
    timer.start();

    vector<std::thread> exec_threads(threadNum);
    int workload = qNum / threadNum;

    for (int thread_index = 0; thread_index < threadNum; thread_index++) {

        exec_threads[thread_index] = std::thread(std::bind([&](const int thread_index) {

            int start = workload * thread_index;
            int end = workload * thread_index + workload;
            end = end > qNum ? qNum: end;

            for (int qID = start; qID < end; qID++) {

                int heapCount = 0;
                const double *qRow = q + qID * d;
                VectorElement *heap = results + qID * k;

                for (int pID = 0; pID < k; pID++) {
                    const double *pRow = p + pID * d;
                    double value = std::inner_product(qRow, qRow + d, pRow, 0.0);
                    heap_enqueue(value, pID, heap, &heapCount);
                }

                for (int pID = k; pID < pNum; pID++) {
                    const double *pRow = p + pID * d;
                    double value = std::inner_product(qRow, qRow + d, pRow, 0.0);

                    if (value > heap[0].data) {
                        heap_dequeue(heap, &heapCount);
                        heap_enqueue(value, pID, heap, &heapCount);
                    }
                }
            }

        }, thread_index));
    }

    for (int thread_index = 0; thread_index < threadNum; thread_index++) {
        exec_threads[thread_index].join();
    }

    timer.stop();

    cout << "naive parallel solution: " << timer.getElapsedTime() << " secs" << endl;
}


int main(int argc, char **argv) {

    bool compareWithNaive = true;
    string QFilePath = "../../data/MovieLens/q.txt";
    string PFilePath = "../../data/MovieLens/p.txt";
    string outputFilePath = "result.txt";
    int k = 10;
    int s = 1000;
    int QNum, PNum, d;
    double *QData = nullptr;
    double *PData = nullptr;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("compareWithNaive", po::value<bool>(&compareWithNaive)->default_value(true), "compare with naive solution")
            ("k", po::value<int>(&k)->default_value(10), "top k")
            ("s", po::value<int>(&s)->default_value(1000), "number of samples")
            ("q_file", po::value<string>(&QFilePath)->default_value("../../data/MovieLens/q.txt"),
             "file path to Q data file")
            ("p_file", po::value<string>(&PFilePath)->default_value("../../data/MovieLens/p.txt"),
             "file path to P data file")
            ("outputFilePath", po::value<string>(&outputFilePath)->default_value("result.txt"),
             "file path to result file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << endl;
        return 0;
    }

    FileUtil::read(QFilePath, QNum, d, QData);
    FileUtil::read(PFilePath, PNum, d, PData);

    cout << "compareWithNaive: " << compareWithNaive << endl;
    cout << "s: " << s << endl;
    cout << "QNum: " << QNum << endl;
    cout << "PNum: " << PNum << endl;
    cout << "d: " << d << endl;

    VectorElement *wedgeSampleResults = new VectorElement[QNum * k];

    WedgeSampling wedgeSampling(QNum, PNum, d, QData, PData);
    wedgeSampling.topk(k, s, wedgeSampleResults);

    if (compareWithNaive) {
        VectorElement *naiveResults = new VectorElement[QNum * k];

        int nthreads = std::thread::hardware_concurrency();
        cout << "number of threads: " << nthreads << endl;
        naive(nthreads, QNum, PNum, d, k, QData, PData, naiveResults);

        FileUtil::outputResult(k, d, QNum, QData, PData, wedgeSampleResults, naiveResults, outputFilePath);
        delete[] naiveResults;
    } else {
        FileUtil::outputResult(k, d, QNum, QData, PData, wedgeSampleResults, outputFilePath);
    }

    delete[] QData;
    delete[] PData;
    delete[] wedgeSampleResults;

    return 0;
}