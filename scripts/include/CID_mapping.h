#ifndef CIDMAPPING_H
#define CIDMAPPING_H

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
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include<time.h>
#include "MurmurHash3.h"
#include <unordered_set>
#include <set>
#include <bitset>
#include <mutex>
#include <condition_variable>
#include <queue>

using namespace std;
namespace CIDmapping {
    string CIDmap_processer(string readFile, string indexFold, int thread, int kmerlen, int bucketnum, string outFile,string thresholdFile);

} // namespace CIDmapping

#endif // CIDMAPPING_H