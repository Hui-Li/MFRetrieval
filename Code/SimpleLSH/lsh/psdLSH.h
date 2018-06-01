//////////////////////////////////////////////////////////////////////////////
/// Copyright (C) 2014 Gefu Tang <tanggefu@gmail.com>. All Rights Reserved.
///
/// This file is part of LSHBOX.
///
/// LSHBOX is free software: you can redistribute it and/or modify it under
/// the terms of the GNU General Public License as published by the Free
/// Software Foundation, either version 3 of the License, or(at your option)
/// any later version.
///
/// LSHBOX is distributed in the hope that it will be useful, but WITHOUT
/// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
/// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
/// more details.
///
/// You should have received a copy of the GNU General Public License along
/// with LSHBOX. If not, see <http://www.gnu.org/licenses/>.
///
/// @version 0.1
/// @author Gefu Tang & Zhifeng Xiao
/// @date 2014.6.30
//////////////////////////////////////////////////////////////////////////////


#ifndef PSDLSH_H
#define PSDLSH_H

#include <random>
#include "../util/Base.h"
#include "../structs/Matrix.h"
#include "../util/Scanner.h"

/**
 * Locality-Sensitive Hashing Scheme Based on p-Stable Distributions.
 *
 *
 * For more information on p-stable distribution based LSH, see the following reference.
 *
 *     Mayur Datar , Nicole Immorlica , Piotr Indyk , Vahab S. Mirrokni,
 *     Locality-sensitive hashing scheme based on p-stable distributions,
 *     Proceedings of the twentieth annual symposium on Computational geometry, June
 *     08-11, 2004, Brooklyn, New York, USA.
 */
class psdLSH
{
public:
    struct Parameter
    {
        /// Hash table size
        unsigned hash_table_size;
        /// Number of hash tables
        unsigned num_of_hash_tables;
        /// Dimension of the vector, it can be obtained from the instance of Matrix
        unsigned D;
        /// Index mode, you can choose 1(CAUCHY) or 2(GAUSSIAN)
        unsigned T;
        /// Window size
        double windows_size;
    };

    psdLSH() {}

    psdLSH(const Parameter &param_)
    {
        init(param_);
    }

    ~psdLSH() {}

    void init(const Parameter &param_);

    void hash(Matrix &data);

   /**
     * Insert a vector to the index.
     *
     * @param key   The sequence number of vector
     * @param ptr The pointer to the vector
     */
    void insert(unsigned key, const double *ptr);

   /**
     * Query the approximate nearest neighborholds.
     *
     * @param ptr   The pointer to the vector
     * @param scanner Top-K scanner, use for scan the approximate nearest neighborholds
     */
    void query(const double *ptr, Scanner<Matrix::Accessor> &scanner);

    /**
     * get the hash value of a vector.
     *
     * @param k     The idx of the table
     * @param ptr The pointer to the vector
     * @return      The hash value
     */
    unsigned getHashVal(unsigned k, const double *ptr);

private:
    Parameter param;
    std::vector<double> rndBs;
    std::vector<std::vector<double> > stableArray;
    std::vector<std::map<unsigned, std::vector<unsigned> > > tables;
};


void psdLSH::init(const Parameter &param_) {
    param = param_;
    tables.resize(param.num_of_hash_tables);
    stableArray.resize(param.num_of_hash_tables);
    std::mt19937 rng(unsigned(std::time(0)));
    std::uniform_real_distribution<double> ur(0, param.windows_size);
    switch (param.T) {
        case CAUCHY: {
            std::cauchy_distribution<double> cd;
            for (std::vector<std::vector<double> >::iterator iter = stableArray.begin();
                 iter != stableArray.end(); ++iter) {
                for (unsigned i = 0; i != param.D; ++i) {
                    iter->push_back(cd(rng));
                }
                rndBs.push_back(ur(rng));
            }
            return;
        }
        case GAUSSIAN: {
            std::normal_distribution<double> nd;
            for (std::vector<std::vector<double> >::iterator iter = stableArray.begin();
                 iter != stableArray.end(); ++iter) {
                for (unsigned i = 0; i != param.D; ++i) {
                    iter->push_back(nd(rng));
                }
                rndBs.push_back(ur(rng));
            }
            return;
        }
        default: {
            return;
        }
    }
}

void psdLSH::hash(Matrix &data) {
    for (unsigned i = 0; i != data.rowNum; ++i) {
        insert(i, data[i]);
    }
}

void psdLSH::insert(unsigned key, const double *ptr) {
    for (unsigned k = 0; k != param.num_of_hash_tables; ++k) {
        unsigned hashVal = getHashVal(k, ptr);
        tables[k][hashVal].push_back(key);
    }
}


void psdLSH::query(const double *ptr, Scanner<Matrix::Accessor> &scanner) {
    scanner.reset(ptr);
    for (unsigned k = 0; k != param.num_of_hash_tables; ++k) {
        unsigned hashVal = getHashVal(k, ptr);
        if (tables[k].find(hashVal) != tables[k].end()) {
            for (std::vector<unsigned>::iterator iter = tables[k][hashVal].begin();
                 iter != tables[k][hashVal].end(); ++iter) {
                scanner(*iter);
            }
        }
    }
    scanner.topk().genTopk();
}

unsigned psdLSH::getHashVal(unsigned k, const double *ptr) {
    double sum(0);
    for (unsigned i = 0; i != param.D; ++i) {
        sum += ptr[i] * stableArray[k][i];
    }
    unsigned hashVal = unsigned(std::floor((sum + rndBs[k]) / param.windows_size)) % param.hash_table_size;
    return hashVal;
}

#endif //PSDLSH_H