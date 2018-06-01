#ifndef BASEUTIL_H
#define BASEUTIL_H

#include <iostream>
#include <fstream>
#include <numeric>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <list>
#include <thread>
#include <string>

using std::cout;
using std::cerr;
using std::endl;
using std::stringstream;
using std::fstream;
using std::ifstream;
using std::ofstream;

using std::string;
using std::vector;
using std::set;
using std::list;
using std::make_pair;
using std::pair;
using std::map;
using std::unordered_map;
using std::unordered_set;
using std::to_string;

////////////////// Boost //////////////////////

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

////////////////// Boost //////////////////////

#define CAUCHY   1
#define GAUSSIAN 2

#define L1_DIST 1
#define L2_DIST 2

#endif //BASEUTIL_H
