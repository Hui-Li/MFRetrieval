#ifndef ARG_MAX_H
#define ARG_MAX_H

#include "../util/Base.h"

struct arg_max {
    double value;
    size_t idx;
    arg_max(size_t idx=0, double value=-1e300): value(value), idx(idx){}
    bool operator<(const arg_max &other) const {return value < other.value;}

    static arg_max default_value() {return arg_max(0, -1e300);}
    static const arg_max& op(const arg_max& a, const arg_max& b) { return a.value < b.value ? b: a; }
};

#endif //ARG_MAX_H