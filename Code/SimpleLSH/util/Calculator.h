#ifndef CALCULATOR_H
#define CALCULATOR_H

#include "../util/Base.h"

inline double sqr(const double x) {
    return x * x;
}

inline double cal_len(const double *vec, int colNum) {
    double len = 0;
    for (int j = 0; j < colNum; ++j) {
        len += vec[j] * vec[j];
    }
    return sqrt(len);
}

#endif //CALCULATOR_H
