#include "util/Base.h"
#include "alg/WedgeNNSampling.h"
#include "util/EvalUtil.h"

/**
 * Wedge Sampling with nonnegativity assumption
 * Diamond Sampling for Approximate Maximum All-pairs Dot-product (MAD) Search, ICDM'15
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {

    string QFilePath = "../../data/MovieLens/q.txt";
    string PFilePath = "../../data/MovieLens/p.txt";
    string groundTruthFilePath = "../../data/MovieLens/result.txt";
    string outputFilePath = "result.txt";
    int k = 10;
    int s = 100;
    bool verify = true;
    bool optimize = true;
    int QNum, PNum, d;
    double *QData = nullptr;
    double *PData = nullptr;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("verify", po::value<bool>(&verify)->default_value(true), "do verification")
            ("optimize", po::value<bool>(&optimize)->default_value(true), "optimized implementation")
            ("k", po::value<int>(&k)->default_value(10), "top k")
            ("s", po::value<int>(&s)->default_value(1000), "number of samples")
            ("q_file", po::value<string>(&QFilePath)->default_value("../../data/MovieLens/q.txt"),
             "file path to Q data file")
            ("p_file", po::value<string>(&PFilePath)->default_value("../../data/MovieLens/p.txt"),
             "file path to P data file")
            ("groundTruth", po::value<string>(&groundTruthFilePath)->default_value("../../data/MovieLens/result.txt"),
             "file path to ground truth")
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

    cout << "k: " << k << endl;
    cout << "s: " << s << endl;
    cout << "optimize: " << optimize << endl;
    cout << "verify: " << verify << endl;
    cout << "QNum: " << QNum << endl;
    cout << "PNum: " << PNum << endl;
    cout << "d: " << d << endl;

    VectorElement *wedgeNNSampleResults = new VectorElement[QNum * k];

    WedgeNNSampling wedgeNNSampling(QNum, PNum, d, QData, PData);
    wedgeNNSampling.topk(k, s, verify, optimize, wedgeNNSampleResults);

    FileUtil::outputResult(k, d, QNum, QData, PData, wedgeNNSampleResults, outputFilePath);

    if (groundTruthFilePath.length() != 0) {
        EvalUtil::avg_recall(QNum, k, groundTruthFilePath, wedgeNNSampleResults);
    }

    delete[] QData;
    delete[] PData;
    delete[] wedgeNNSampleResults;

    return 0;
}