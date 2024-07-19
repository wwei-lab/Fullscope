#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <ranges>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <future>
#include <thread>
#include <map>
#include <sstream>
#include <inttypes.h>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/search/search.hpp>
#include <algorithm> // For std::transform (optional, but can be used)  
#include <cctype> 

using namespace std;
using namespace seqan3;
using namespace seqan3::literals;
vector<string> anchor{ "ATGGCGACCTTATCAG", "TTGTCTTCCTAAGACCG" };
auto left_dna5 = anchor[0] | seqan3::views::char_to<seqan3::dna5>;
auto right_dna5 = anchor[1] | seqan3::views::char_to<seqan3::dna5>;
auto method = seqan3::align_cfg::method_local{};
seqan3::align_cfg::scoring_scheme scheme{ seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}} };
seqan3::align_cfg::gap_cost_affine gap_costs{ seqan3::align_cfg::open_score{0}, seqan3::align_cfg::extension_score{-1} };
auto aligncfg = method | scheme | gap_costs;

vector<int> mywidth{ 20, 30 , 25 };
struct ReadInfo {
    string id;
    string seq;
    string gene;
};
std::string reverseComplement(std::string mStr) {
    std::string str(mStr.length(), 0);
    for (int c = 0; c < mStr.length(); c++) {
        char base = mStr[c];
        switch (base) {
        case 'A':
        case 'a':
            str[mStr.length() - c - 1] = 'T';
            break;
        case 'T':
        case 't':
            str[mStr.length() - c - 1] = 'A';
            break;
        case 'C':
        case 'c':
            str[mStr.length() - c - 1] = 'G';
            break;
        case 'G':
        case 'g':
            str[mStr.length() - c - 1] = 'C';
            break;
        default:
            str[mStr.length() - c - 1] = 'N';
        }
    }
    return str;
}

void loadGeneRanges(string gtfFile, vector<pair <string, string>>& geneRanges) {
    ifstream inFile(gtfFile);
    if (!inFile) {
        cerr << "Error: failed to open " << gtfFile << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    while (getline(inFile, line)) {
        istringstream ss(line);
        string seqname, genename;
        string startc, endc;
        uint64_t start, end;
        getline(ss, seqname, '\t');
        getline(ss, genename, '\t');
        getline(ss, startc, '\t');
        getline(ss, endc, '\t');
        start = stoull(startc);
        end = stoull(endc);

        if (start > 1000) {
            start = start - 1000;
        }
        else {
            start = 0;
        }

        end = end + 1000;
        startc = to_string(start);
        endc = to_string(end);
        string rg = seqname + ":" + startc + "-" + endc;
        geneRanges.push_back(make_pair(genename, rg));
    }
}


struct HitInfo {
    int forward1 = -1;
    int forward2 = -1;
    int reverse1 = -1;
    int reverse2 = -1;
};

std::string dna5_vector_substring_to_string(const seqan3::dna5_vector& vec, size_t start, size_t end) {
    // 确保索引在有效范围内  
    if (start >= vec.size() || end > vec.size() || start > end) {
        throw std::out_of_range("Indices out of range");
    }

    std::string result;
    result.reserve(end - start); // 可选，预先分配内存  
    for (size_t i = start; i < end; ++i) {
        result += seqan3::to_char(vec[i]); // 使用 seqan3 提供的函数将 dna5 转换为 char  
    }
    return result;
}


seqan3::dna5_vector string_to_dna5_vector(const std::string& s) {
    seqan3::dna5_vector result;
    result.reserve(s.size()); // 可选，预先分配内存以提高效率  
    for (char c : s) {
        if (c == 'A') result.push_back('A'_dna5);
        else if (c == 'C') result.push_back('C'_dna5);
        else if (c == 'G') result.push_back('G'_dna5);
        else if (c == 'T') result.push_back('T'_dna5);
        else if (c == 'N') result.push_back('N'_dna5);
    }
    return result;
}

string extractCIDs(std::vector<dna5_vector>& seqvec, std::vector<string>& readids,string geneid,
    std::map<uint32_t, HitInfo>& hitinfo_map){
    //std::vector<std::string> cids(seqvec.size(), ""); // 存储提取的 CID 序列
    std::string cid = "";
    std::string strand = "";
    string outlines = "";
    // 遍历每一个序列
    for (size_t i = 0; i < seqvec.size(); ++i) {
        auto hitinfo = hitinfo_map[i];
        bool has_forward_hit = (hitinfo.forward1 >= 0 || hitinfo.forward2 >= 0);
        bool has_reverse_hit = (hitinfo.reverse1 >= 0 || hitinfo.reverse2 >= 0);

        if (has_forward_hit && !has_reverse_hit) {
            // 仅有正向匹配，将其往相应方向延伸25bp作为 reference id
            uint64_t start_pos = (hitinfo.forward1 >= 0) ? hitinfo.forward1 + anchor[0].size(): hitinfo.forward2 - 25;
            if (start_pos < 0) start_pos = 0;
            uint64_t end_pos = (hitinfo.forward2 >= 0) ? hitinfo.forward2 : hitinfo.forward1 + 25;
            if (end_pos > seqvec[i].size()) end_pos = seqvec[i].size();
            if (end_pos > start_pos) {
                cid = dna5_vector_substring_to_string(seqvec[i], start_pos, end_pos);
                strand = "+";
            }
            //cid = std::string{ seqvec[i].begin() + start_pos, seqvec[i].begin() + end_pos };
        }
        else if (!has_forward_hit && has_reverse_hit) {
            // 仅有反向匹配，将其反向互补并往相应方向延伸25bp作为 reference id
            uint64_t start_pos = (hitinfo.reverse2 >= 0) ? hitinfo.reverse2 + anchor[1].size(): hitinfo.reverse1 - 25;
            if (start_pos < 0) start_pos = 0;
            uint64_t end_pos = (hitinfo.reverse1 >= 0) ? hitinfo.reverse1 : hitinfo.reverse2 + 25;
            if (end_pos > seqvec[i].size()) end_pos = seqvec[i].size();
            if (end_pos > start_pos) {
                cid = dna5_vector_substring_to_string(seqvec[i], start_pos, end_pos);
                cid = reverseComplement(cid);
                strand = "-";
            }
        }
        else if (has_forward_hit && has_reverse_hit) {
            // 如果正向和反向匹配都存在，则舍弃该序列
            cid = "discard";
        }
        else {
            // 如果都没有匹配，则返回 none
            cid = "none";
        }
        outlines += readids[i] + "\t*\t" + cid + "\t" + geneid + "\t" + strand + "\n";
    }
    return outlines;
}

/*
* std::string reverseComplement(dna5& seq) {
    auto seq_dna5 = seq | views::char_to<dna5>;
    auto rc_seq_dna5 = seq_dna5 | views::complement | views::reverse;
    std::string rc_seq = rc_seq_dna5 | views::to<std::string>;
    return rc_seq;
}
*/


void handleBamChunk(string bamFile, uint64_t gstart, uint64_t gend, vector<pair <string, string>>& geneRanges,
    uint16_t errorcount, ofstream& outputFile, mutex& mtx) {
    //cout << "open bam file" << endl;
    samFile* fp = sam_open(bamFile.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(fp);
    hts_idx_t* idx = sam_index_load(fp, bamFile.c_str());
    bam1_t* record = bam_init1();
    auto search_cfg = seqan3::search_cfg::max_error_total{ seqan3::search_cfg::error_count{(unsigned char)errorcount} } | seqan3::search_cfg::hit_all_best{};

    auto left_dna5_rc = reverseComplement(anchor[0]) | seqan3::views::char_to<seqan3::dna5>;
    auto right_dna5_rc = reverseComplement(anchor[1]) | seqan3::views::char_to<seqan3::dna5>;

    for (int gi = gstart; gi <= gend; gi++) {
        string geneid = geneRanges[gi].first;
        string rg = geneRanges[gi].second;
        const char* c_rg = rg.c_str();
        hts_itr_t* iter = sam_itr_querys(idx, header, c_rg);
        std::vector<dna5_vector> seqvec;
        std::vector<string> readids;

        while (sam_itr_next(fp, iter, record) >= 0) {
            string readid = bam_get_qname(record);
            uint8_t* seq = bam_get_seq(record);
            std::string seq_str;

            for (int i = 0; i < record->core.l_qseq; ++i) {
                seq_str.push_back(seq_nt16_str[bam_seqi(seq, i)]);
                //seqan3::dna5 d2 = seqan3::assign_char_to(chari, seqan3::dna5{});
                //seq_dna5.push_back(d2);
            }
           // seqan3::dna5_vector seq_dna5 = ;
            uint64_t pos = record->core.pos;
            uint64_t end = bam_endpos(record);
            seqvec.push_back(string_to_dna5_vector(seq_str));
            readids.push_back(readid);
        }
        if (seqvec.size() == 0) { continue; }
        hts_itr_destroy(iter);
        seqan3::bi_fm_index indexi{ seqvec };
        auto results1 = search(left_dna5, indexi, search_cfg);
        auto results2 = search(right_dna5, indexi, search_cfg);
        auto results_rc1 = search(left_dna5_rc, indexi, search_cfg);
        auto results_rc2 = search(right_dna5_rc, indexi, search_cfg);

        std::map<uint32_t, HitInfo> hitinfo_map;
        for (auto& res : results1) {
            uint32_t refid = res.reference_id();
            uint32_t hitpos = res.reference_begin_position();
            hitinfo_map[refid].forward1 = (hitpos);
        }
        for (auto& res : results2) {
            uint32_t refid = res.reference_id();
            uint32_t hitpos = res.reference_begin_position();
            hitinfo_map[refid].forward2 = (hitpos);
        }
        for (auto& res : results_rc1) {
            uint32_t refid = res.reference_id();
            uint32_t hitpos = res.reference_begin_position();
            hitinfo_map[refid].reverse1 = (hitpos);
        }
        for (auto& res : results_rc2) {
            uint32_t refid = res.reference_id();
            uint32_t hitpos = res.reference_begin_position();
            hitinfo_map[refid].reverse2 = (hitpos);
        }
        //bi_fm_index<dna5>& index;
        string outlines = extractCIDs(seqvec, readids, geneid, hitinfo_map);
        mtx.lock();
        outputFile << outlines;
        mtx.unlock();
    }

    bam_destroy1(record);
    hts_idx_destroy(idx);
    sam_close(fp);
    bam_hdr_destroy(header);

}

void handleBam(string bamFile, const string& outputFileName, vector<pair <string, string>>& geneRanges, uint16_t errorcount, int numThreads) {
    //cout << "bam open " << endl;
    ofstream outputFile(outputFileName);
    vector<thread> threads;
    mutex mtx;
    uint32_t genecount = geneRanges.size();
    uint32_t chunk_size = genecount / numThreads;
    uint32_t remainder = genecount % numThreads;
    uint32_t start = 0, end = 0;
    for (uint32_t i = 0; i < numThreads; i++) {
        start = end;
        end = start + chunk_size;
        if (i == numThreads - 1) {
            end += remainder;
        }
        cout << start << "-" << end << endl;
        threads.emplace_back(handleBamChunk, bamFile, start, end, ref(geneRanges), errorcount, ref(outputFile), ref(mtx));
    }

    for (auto& t : threads) {
        t.join();
    }
}


int main(int argc, char* argv[]) {
    cout << "usage:" << endl;
    string bamFileA = argv[1];
    string geneAn = argv[2];
    string outFile = argv[3];
    int thread = atoi(argv[4]);
    uint16_t errorcount = atoi(argv[5]);

    //double errorate = atof(argv[4]);
    //part of cidprimer, capture oligo
    cout << "load bam: " << endl;
    //vector<ReadInfo> reads;
    vector<pair <string, string>> geneRanges;
    loadGeneRanges(geneAn, geneRanges);
    handleBam(bamFileA, outFile, geneRanges, errorcount, thread);
    //findCIDSequenceThread(reads,  outprex, thread, errorcount);

    return 0;
}