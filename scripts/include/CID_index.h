// CIDindex.h
#ifndef CIDINDEX_H
#define CIDINDEX_H

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include "MurmurHash3.h"
#include <time.h>
#include <map>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/bitset.hpp>  // for std::bitset
#include <htslib/sam.h>
#include <set>
#include <thread>
#include <mutex>
#include <sstream>
#include <unistd.h> 
#include <bitset>
using namespace std;
namespace CIDindex {
    // 配置参数建议通过结构体传递
    void build_index(string tableFile, int kmerlen, int bucketnum, int numThreads, string outFold);
}

#endif