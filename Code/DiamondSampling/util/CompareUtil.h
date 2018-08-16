#ifndef COMPAREUTIL_H
#define COMPAREUTIL_H

#include "Base.h"

namespace CompareUtil {
    // This function returns true if the first pair is "larger"
    // than the second one according to some metric
    template<typename T>
    bool pairGreaterCompare(const std::pair<T, double> &firstElem, const std::pair<T, double> &secondElem) {
        return firstElem.second > secondElem.second;
    }

    // This function returns true if the first pair is "less"
    // than the second one according to some metric
    template<typename T>
    bool pairLessCompare(const std::pair<T, double> &firstElem, const std::pair<T, double> &secondElem) {
        return firstElem.second < secondElem.second;
    }

    // descending comparator for arg_sort
    template<typename T>
    struct comparator {
        const T* value;
        comparator(const T* value=NULL): value(value){}
        bool operator()(size_t a, size_t b) const { return value[a] > value[b]; }
    };
}

#endif //COMPAREUTIL_H
