#include "util/Base.h"
#include "util/Conf.h"
#include "structs/Matrix.h"
#include "util/Logger.h"
#include "alg/Naive.h"
#include "alg/tree/BallTreeSearch.h"
#include "alg/svd/SVDIncrPrune.h"
#include "alg/svd/SVDIntUpperBoundIncrPrune.h"
#include "alg/svd/SVDIntUpperBoundIncrPrune2.h"
#include "alg/int/IntUpperBound.h"
#include "alg/tree/FastMKS.h"
#include "alg/transformation/TransformIncrPrune.h"
#include "alg/simd/SIMDIntUpperBound.h"
#include "alg/svd/SVDIncrPruneIndividualReorder.h"
#include "alg/transformation/TransformSVDIncrPrune.h"
#include "alg/transformation/TransformSVDIncrPrune2.h"
#include "alg/int/IntUpperBound2.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

void basicLog(const Matrix &q, const Matrix &p, const int k) {
    Logger::Log("q path: " + to_string(Conf::qDataPath));
    Logger::Log("p path: " + to_string(Conf::pDataPath));
    Logger::Log("q: " + to_string(q.rowNum) + "," + to_string(q.colNum));
    Logger::Log("p: " + to_string(p.rowNum) + "," + to_string(p.colNum));
    Logger::Log("Algorithm: " + Conf::algName);
    Logger::Log("k: " + to_string(k));
}

int main(int argc, char **argv) {

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("alg", po::value<string>(&(Conf::algName))->default_value("naive"), "Algorithm")
            ("k", po::value<int>(&(Conf::k))->default_value(1), "K")
            ("dataset", po::value<string>(&(Conf::dataset)), "name of dataset for log output")
            ("q", po::value<string>(&(Conf::qDataPath)), "file path of q Data")
            ("p", po::value<string>(&(Conf::pDataPath)), "file path of p Data")
            ("scalingValue", po::value<int>(&(Conf::scalingValue))->default_value(127), "maximum value for scaling")
            ("SIGMA", po::value<double>(&(Conf::SIGMA))->default_value(0.8), "percentage of SIGMA for SVD incremental prune")
            ("log", po::value<bool>(&(Conf::log))->default_value(true), "whether it outputs log")
            ("logPathPrefix", po::value<string>(&(Conf::logPathPrefix))->default_value("./log"), "output path of log file (Prefix)")
            ("outputResult", po::value<bool>(&(Conf::outputResult))->default_value(true), "whether it outputs results")
            ("resultPathPrefix", po::value<string>(&(Conf::resultPathPrefix))->default_value("./result"), "output path of result file (Prefix)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 0;
    } else if (Conf::qDataPath == "" || Conf::pDataPath == "") {
        cout << "Please specify path to data files" << endl << endl;
        cout << desc << endl;
        return 0;
    }

//    Conf::qDataPath = "../../data/MovieLens/q.txt";
//    Conf::pDataPath = "../../data/MovieLens/p.txt";
//    Conf::dataset = "MovieLens";
//    Conf::k = 1;
//    Conf::SIGMA = 0.7;
//    Conf::algName = "FEIPR-I2";
//    Conf::algName = "FEIPR-I";
//    Conf::algName = "FEIPR-S";

    Conf::Output();

    Matrix q;
    Matrix p;
    q.readData(Conf::qDataPath);
    p.readData(Conf::pDataPath);
    Conf::dimension = p.colNum;

    cout << "-----------------------" << endl;

    // ToDo: replace the old name (FEIPR) with FEXIPRO

    if (Conf::algName == "Naive") {
        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + ".txt";
        Logger::open(logFileName);
        basicLog(q, p, Conf::k);

        naive(Conf::k, q, p);

    } else if (Conf::algName == "SIR") {

        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + "-" + to_string(Conf::scalingValue) + "-" + to_string(Conf::SIGMA) + ".txt";
        Logger::open(logFileName);
        basicLog(q, p, Conf::k);
        Logger::Log("SIGMA: " + to_string(Conf::SIGMA));
        Logger::Log("Scaling Value: " + to_string(Conf::scalingValue));

        SIRPrune sirPrune(Conf::k, Conf::scalingValue, Conf::SIGMA, &q, &p);
        sirPrune.topK();

    } else if (Conf::algName == "SR") {

        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + "-" + to_string(Conf::scalingValue) + "-" + to_string(Conf::SIGMA) + ".txt";
        Logger::open(logFileName);
        basicLog(q, p, Conf::k);
        Logger::Log("SIGMA: " + to_string(Conf::SIGMA));

        TransformSVDIncrPrune2 transformSVDIncrPrune2(Conf::k, Conf::SIGMA, &q, &p);
        transformSVDIncrPrune2.topK();

    } else if (Conf::algName == "I") {

        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + "-" + to_string(Conf::scalingValue) + ".txt";
        Logger::open(logFileName);
        basicLog(q, p, Conf::k);
        Logger::Log("Scaling Value: " + to_string(Conf::scalingValue));

        IntUpperBound intUpperBound(Conf::k, Conf::scalingValue, &q, &p);
        intUpperBound.topK();

    } else if (Conf::algName == "I-SIMD") {

        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + "-" + to_string(Conf::scalingValue) + ".txt";
        Logger::open(logFileName);
        basicLog(q, p, Conf::k);
        Logger::Log("Scaling Value: " + to_string(Conf::scalingValue));

        SIMDIntUpperBound simdIntUpperBound(Conf::k, &q, &p);
        simdIntUpperBound.topK();

    } else if (Conf::algName == "S"){
        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset  + "-" + Conf::algName + "-" + to_string(Conf::k) + "-" +
                             to_string(Conf::SIGMA) + ".txt";
        Logger::open(logFileName);
        basicLog(q, p, Conf::k);
        Logger::Log("SIGMA: " + to_string(Conf::SIGMA));

        SVDIncrPrune svdIncrPrune(Conf::k, Conf::SIGMA, &q, &p);
        svdIncrPrune.topK();

    } else if (Conf::algName == "S-Ind"){
        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset  + "-" + Conf::algName + "-" + to_string(Conf::k) + "-" +
                             to_string(Conf::SIGMA) + ".txt";
        Logger::open(logFileName);
        basicLog(q, p, Conf::k);
        Logger::Log("SIGMA: " + to_string(Conf::SIGMA));

        SVDIncrPruneIndividualReorder svdIncrPruneIndividualReorder(Conf::k, Conf::SIGMA, &q, &p);
        svdIncrPruneIndividualReorder.topK();

    } else if (Conf::algName == "SI") {

        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset  + "-" + Conf::algName + "-" + to_string(Conf::k) + "-" + to_string(Conf::scalingValue) + "-" +
                             to_string(Conf::SIGMA) + ".txt";
        Logger::open(logFileName);
        basicLog(q, p, Conf::k);
        Logger::Log("SIGMA: " + to_string(Conf::SIGMA));
        Logger::Log("Scaling Value: " + to_string(Conf::scalingValue));

        SVDIntUpperBoundIncrPrune svdIntUpperBoundIncrPrune(Conf::k, Conf::scalingValue, Conf::SIGMA, &q, &p);
        svdIntUpperBoundIncrPrune.topK();

    } else if (Conf::algName == "BallTree") {

        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset  + "-" + Conf::algName + "-" + to_string(Conf::k) + ".txt";
        Logger::open(logFileName);
        basicLog(q, p, Conf::k);

        ballTreeTopK(Conf::k, q, p);

    } else if (Conf::algName == "FastMKS") {

        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset  + "-" + Conf::algName + "-" + to_string(Conf::k) + ".txt";
        Logger::open(logFileName);
        basicLog(q, p, Conf::k);

        fastMKS(Conf::k, Conf::pDataPath, Conf::qDataPath);
    } else {
        cout << "unrecognized method" << endl;
    }

    return 0;
}
