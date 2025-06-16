#pragma once
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <bitset>

struct AlignmentResult {
    int start;
    int end;
    int distance;
    std::string adapter_id;
    bool is_rc;
};

struct Segment {
    int start;
    int end;
    std::string adapter;
    std::string strand;
    std::string read_id;
};

struct ReadData {
    std::string header;
    std::string sequence;
    std::string quality;
    size_t read_num;
};

struct GeneRange {
    std::string chr;
    uint64_t start;
    uint64_t end;
};

struct SeqInfo {
    std::bitset<50> encode_seqs;  // K=25, 2 bits per base
    uint32_t xpos;
    uint32_t ypos;
};