#include "util/Base.h"
#include "structs/VectorElement.h"
#include "util/Monitor.h"
#include "structs/FastHeap.h"
#include "util/FileUtil.h"
#include "alg/BudgetedMIP.h"
#include "util/EvalUtil.h"


int main(int argc, char **argv) {
    string QFilePath = "../../data/MovieLens/q.txt";
    string PFilePath = "../../data/MovieLens/p.txt";
    string groundTruthFilePath = "../../data/MovieLens/result.txt";
    string outputFilePath = "naive_result.txt";
    int k = 10;
    int budget = 1024;
    int QNum, PNum, d;
    double *QData = nullptr;
    double *PData = nullptr;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("k", po::value<int>(&k)->default_value(10), "top k")
            ("b", po::value<int>(&budget)->default_value(1024), "budget")
            ("q_file", po::value<string>(&QFilePath)->default_value("../../data/MovieLens/q.txt"),
             "file path to Q data file")
            ("p_file", po::value<string>(&PFilePath)->default_value("../../data/MovieLens/p.txt"),
             "file path to P data file")
            ("groundTruth", po::value<string>(&groundTruthFilePath)->default_value("../../data/MovieLens/result.txt"),
             "file path to ground truth")
            ("outputFilePath", po::value<string>(&outputFilePath)->default_value("naive_result.txt"),
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

    cout << "k: " << k << endl;
    cout << "Budget: " << budget << endl;
    cout << "QNum: " << QNum << endl;
    cout << "PNum: " << PNum << endl;
    cout << "d: " << d << endl;

    Monitor timer;

    timer.start();
    vector<vector<size_t> > results(QNum);
    BudgetedMIP budgetedMIP(QNum, PNum, d, QData, PData);
    budgetedMIP.topk(k, budget, results);
    timer.stop();

    cout << "retrieval time: " << timer.getElapsedTime() << endl;

    int nthreads = std::thread::hardware_concurrency();

    EvalUtil::avg_recall(QNum, k, groundTruthFilePath, results);
    FileUtil::outputResult(k, d, QNum, QData, PData, results, outputFilePath);

    delete[] QData;
    delete[] PData;

    return 0;
}