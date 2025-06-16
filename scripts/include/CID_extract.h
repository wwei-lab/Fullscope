// cidmap_extract.h
#ifndef CIDMAP_EXTRACT_H
#define CIDMAP_EXTRACT_H
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
namespace CIDExtract {

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

using GeneRanges = std::map<std::string, 
    std::map<std::string, std::pair<uint64_t, uint64_t>>>;

// Main interface
string process_bam(const std::string& bam_path,
                const std::string& output_path,
                GeneRanges& gene_ranges,
                Config& config);

// Utility functions
void load_gene_ranges(const std::string& gtf_path, GeneRanges& gene_ranges);
uint64_t count_bam_records(const std::string& bam_path);

} // namespace CIDExtract

#endif // CIDMAP_EXTRACT_H