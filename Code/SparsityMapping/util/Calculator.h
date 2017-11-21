#ifndef CALCULATOR_H
#define CALCULATOR_H

#include "../util/Base.h"

inline int sign(const double x){
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

inline double cal_len(const double *vec, int colNum) {
    double len = 0;
    for (int j = 0; j < colNum; ++j) {
        len += vec[j] * vec[j];
    }
    return sqrt(len);
}

inline double inner_product(const double *a, const double *b, int colNum){
    double result = 0.0;
    for(int i=0;i<colNum;i++){
        result += a[i] * b[i];
    }
    return result;
}

inline double inner_product(double *a, double *b, int colNum){
    double result = 0.0;
    for(int i=0;i<colNum;i++){
        result += a[i] * b[i];
    }
    return result;
}

#endif //CALCULATOR_H
