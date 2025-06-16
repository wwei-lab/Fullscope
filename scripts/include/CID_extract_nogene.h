// cidmap_extract.h
#ifndef CIDMAP_EXTRACT_NOGENE_H
#define CIDMAP_EXTRACT_NOGENE_H
#include <htslib/sam.h>
#include <seqan3/core/all.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alignment/all.hpp>
#include <string>
#include <vector>
#include <map>
#include <mutex>
#include <utility>
#include <thread>
#include <atomic>
using namespace std;
namespace CIDExtractFastq {

struct Config {
    std::vector<std::string> anchor_seqs {
        //"ATGGCGACCTTATCAG", 
        //"TTGTCTTCCTAAGACCG"
        "CGGTCTTAGGAAGACAA",
        "CTGATAAGGTCGCCAT"
        
    };
    
    std::vector<int> width_params {20, 30, 25};
    int error_threshold = 5;
    int num_threads = 4;
};

struct ReadInfo {
    std::string id;
    std::string seq;
    std::string gene;
};

// Main interface
string extractcid_fastq(const string& fastq_path,
                   const string& output_path,
                   Config& config);

} // namespace CIDExtract

#endif // CIDMAP_EXTRACT_H