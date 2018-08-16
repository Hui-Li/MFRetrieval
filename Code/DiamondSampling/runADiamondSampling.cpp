#include "util/Base.h"
#include "structs/VectorElement.h"
#include "util/FileUtil.h"
#include "util/EvalUtil.h"
#include "alg/DiamondSampling.h"
#include "alg/ADiamondSampling.h"

int main(int argc, char **argv) {

    string QFilePath = "../../data/MovieLens/q.txt";
    string PFilePath = "../../data/MovieLens/p.txt";
    string groundTruthFilePath = "../../data/MovieLens/result.txt";
    string outputFilePath = "result.txt";

    int k = 10;
    int s = 1000;
    int QNum, PNum, d;
    double *QData = nullptr;
    double *PData = nullptr;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
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
    cout << "QNum: " << QNum << endl;
    cout << "PNum: " << PNum << endl;
    cout << "d: " << d << endl;

    VectorElement *diamondSampleResults = new VectorElement[QNum * k];

    Monitor timer;
    timer.start();

    ADiamondSampling adiamondSampling(QNum, PNum, d, QData, PData);
    adiamondSampling.topk(k, s, diamondSampleResults);

    timer.stop();

    cout << "time: " << timer.getElapsedTime() << " seconds" << endl;

    FileUtil::outputResult(k, d, QNum, QData, PData, diamondSampleResults, outputFilePath);

    if (groundTruthFilePath.length() != 0) {
        avg_recall(QNum, k, groundTruthFilePath, diamondSampleResults);
    }

    delete[] QData;
    delete[] PData;
    delete[] diamondSampleResults;

    return 0;
}