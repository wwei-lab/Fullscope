#ifndef CIDMAPPINGPRECISE_H
#define CIDMAPPINGPRECISE_H

#include "utils.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <ranges>
#include <future>
#include <thread>
#include <map>
#include <sstream>
#include <cstdio>
#include "edlib.h"
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/bitset.hpp>  // for std::bitset
#include<time.h>
#include "MurmurHash3.h"
#include <unordered_set>
#include <set>
#include <bitset>
#include <mutex>
#include <condition_variable>
#include <queue>

using namespace std;
namespace CIDmappingPrecise {
    string CIDmap_processer_precise(string readFile, string indexFold, int thread, int kmerlen, string outFile);
} // namespace CIDmapping

#endif // CIDMAPPING_H