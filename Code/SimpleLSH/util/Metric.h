#ifndef METRIC_H
#define METRIC_H

#include "Calculator.h"

/**
 * Use for common distance functions.
 */

class Metric {
    unsigned type_;
    unsigned dim_;
public:
    /**
     * Constructor for this class.
     *
     * @param dim  Dimension of each vector
     * @param type The way to measure the distance, you can choose 1(L1_DIST) or 2(L2_DIST)
     */
    Metric(unsigned dim, unsigned type) : dim_(dim), type_(type) {}

    ~Metric() {}

    /**
     * Get the dimension of the vectors
     */
    unsigned dim() const {
        return dim_;
    }

    /**
     * measure the distance.
     *
     * @param  vec1 The first vector
     * @param  vec2 The second vector
     * @return      The distance
     */
    double dist(const double *vec1, const double *vec2) const {
        float dist_ = 0.0;
        switch (type_) {
            case L1_DIST: {
                for (unsigned i = 0; i != dim_; ++i) {
                    dist_ += double(std::abs(vec1[i] * 1.0 - vec2[i]));
                }
                return dist_;
            }
            case L2_DIST: {
                for (unsigned i = 0; i != dim_; ++i) {
                    dist_ += sqr(vec1[i] - vec2[i]);
                }
                return std::sqrt(dist_);
            }
            default: {
                return -1;
            }
        }
    }
};

#endif //METRIC_H
