#ifndef BASE_H
#define BASE_H

#include <iostream>
#include <stdint.h>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>
#include <time.h>
#include <queue>
#include <cmath>
#include <map>
#include <stdexcept>
#include <algorithm>
#include <cfloat>
#include <cstdlib>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

using namespace std;

string getCurrentTime(const char *pattern = "%Y%m%d-%H%M%S") {
    time_t cur_time = time(0);
    char buf[100];
    strftime(buf, 100, pattern, localtime(&cur_time)); //format date and time.
    return buf;
}

template<typename T>
std::string to_string(const T &n) {
    std::ostringstream stm;
    stm << n;
    return stm.str();
}

void split(const string &s, vector<string> &elems, char delim=',') {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

typedef unsigned long long ull;

//#define TIME_IT

#endif //BASE_H
