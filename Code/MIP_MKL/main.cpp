#include "util/Base.h"
#include "util/Conf.h"
#include "util/Logger.h"
#include "alg/ArmaNaive.h"
#include "alg/MKLNaive.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;


int main(int argc, char **argv) {
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("alg", po::value<string>(&(Conf::algName))->default_value("MKLBatchNaive"), "Algorithm")
            ("k", po::value<int>(&(Conf::k))->default_value(1), "K")
            ("batchSize", po::value<int>(&(Conf::batchSize))->default_value(1000), "query number for one batch")
            ("dataset", po::value<string>(&(Conf::dataset)), "name of dataset for log output")
            ("q", po::value<string>(&(Conf::qDataPath))->default_value("../data/q.txt"), "file path of q Data")
            ("p", po::value<string>(&(Conf::pDataPath))->default_value("../data/p.txt"), "file path of p Data")
            ("log", po::value<bool>(&(Conf::log))->default_value(true), "whether output log")
            ("singleThread", po::value<bool>(&(Conf::singleThread))->default_value(false), "single thread")
            ("logPathPrefix", po::value<string>(&(Conf::logPathPrefix))->default_value("./log"), "output path of log file (Prefix)")
            ("outputResult", po::value<bool>(&(Conf::outputResult))->default_value(false), "whether output results")
            ("resultPathPrefix", po::value<string>(&(Conf::resultPathPrefix))->default_value("./result"), "output path of result file (Prefix)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 0;
    }

    Conf::Output();

    Logger::Log("q path: " + to_string(Conf::qDataPath));
    Logger::Log("p path: " + to_string(Conf::pDataPath));
    Logger::Log("Algorithm: " + Conf::algName);
    Logger::Log("k: " + to_string(Conf::k));

    if (Conf::algName == "ArmaNaive") {
        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + ".txt";
        Logger::open(logFileName);

        ArmaNaive::naive();

    } else if (Conf::algName == "ArmaBatchNaive") {

        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + "-" + to_string(Conf::batchSize) + ".txt";
        Logger::open(logFileName);
        Logger::Log("batchSize: " + to_string(Conf::batchSize));

        ArmaNaive::batchNaive();

    } else if (Conf::algName == "MKLNaive") {

        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k);

        if (Conf::singleThread) {
            logFileName += "-singlethread.txt";
        } else {
            logFileName += ".txt";
        }

        Logger::open(logFileName);
        MKLNaive::naive();
    } else if (Conf::algName == "MKLBatchNaive") {

        string logFileName = Conf::logPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + "-" + to_string(Conf::batchSize);

        if (Conf::singleThread) {
            logFileName += "-singlethread.txt";
        } else {
            logFileName += ".txt";
        }

        Logger::open(logFileName);
        Logger::Log("batchSize: " + to_string(Conf::batchSize));

        MKLNaive::batchNaive();
    }

    return 0;
}