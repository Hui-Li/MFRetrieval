#ifndef COMPAREUTIL_H
#define COMPAREUTIL_H

#include "Base.h"

namespace CompareUtil {
    // This function returns true if the first pair is "larger"
    // than the second one according to some metric
    bool pairGreaterCompare(const std::pair<int, double> &firstElem, const std::pair<int, double> &secondElem) {
        return firstElem.second > secondElem.second;
    }

    // This function returns true if the first pair is "less"
    // than the second one according to some metric
    bool pairLessCompare(const std::pair<int, double> &firstElem, const std::pair<int, double> &secondElem) {
        return firstElem.second < secondElem.second;
    }
}

#endif //COMPAREUTIL_H
