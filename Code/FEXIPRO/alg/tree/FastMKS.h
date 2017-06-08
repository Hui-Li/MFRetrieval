#ifndef FASTMKS_H
#define FASTMKS_H

#include "../../util/Monitor.h"
#include <mlpack/methods/fastmks/fastmks_model.hpp>

using namespace mlpack::fastmks;
using namespace mlpack::kernel;
using namespace mlpack::data;
using namespace std;

void fastMKS(const int k, const string pFile, const string qFile) {

    arma::mat referenceData;
    arma::mat queryData;

    bool single = true;
    bool naive = false;
    double base = 1.3;

    Load(qFile, queryData, false, true);
    Load(pFile, referenceData, false, true);

    Logger::Log("--------------------------");
    Logger::Log("FastMKS");
    Logger::Log("Loaded q data from '" + qFile + "' (" + to_string(queryData.n_rows) + " x " + to_string(queryData.n_cols) + ").");
    Logger::Log("Loaded p data from '" + pFile + "' (" + to_string(referenceData.n_rows) + " x " + to_string(referenceData.n_cols) + ").");

    Monitor t;
    t.start();

    arma::Mat<size_t> indices;
    arma::mat products;

    FastMKSModel model;
    model.KernelType() = FastMKSModel::LINEAR_KERNEL;
    LinearKernel lk;
    model.BuildModel(referenceData, lk, single, naive, base);

    t.stop();

    Logger::Log("preprocess: " + to_string(t.getElapsedTime()) + " secs");

    t.start();

    model.Search(queryData, (size_t) k, indices, products, base);

    t.stop();

    Logger::Log("search time: " + to_string(t.getElapsedTime()) + " secs");


    if (Conf::outputResult) {
        string resultFileName = Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + "-result.txt";
        string indexFileName = Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + "-index.txt";
        Save(resultFileName, products, false);
        Save(indexFileName, indices, false);
    }

}

#endif //FASTMKS_H