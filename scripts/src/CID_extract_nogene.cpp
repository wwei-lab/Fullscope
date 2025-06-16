// cidmap_extract.cpp
#include "CID_extract_nogene.h"
#include "utils.h"
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>
#include <memory>
#include <iostream>
#include <string>
#include <unordered_map>
#include <seqan3/io/sequence_file/all.hpp>

using namespace std;
using namespace FsUtils;
namespace CIDExtractFastq {

using namespace seqan3::literals;
uint64_t total_processed = 0;
uint64_t valid_cid_count = 0;

// 保留原有的CidExtractPairwiseAlign函数不变
std::string CidExtractPairwiseAlignKmer(std::string sequence, Config& config, string& align_status) {
    string cid = "";
    const auto& anchors = config.anchor_seqs;
    const auto& widths = config.width_params;
    const int error_count = config.error_threshold;

    // 1. 截取序列末端150bp（或更短）作为目标序列
    const int sub_len = std::min(150, static_cast<int>(sequence.size()));
    const std::string seq = sequence.substr(sequence.size()-sub_len, sub_len);
    const std::string rc_seq = FsUtils::reverse_complement(sequence).substr(sequence.size()-sub_len, sub_len);

    // 2. 建立k-mer索引（k=10）
    const int k = 10;
    auto build_kmer_index = [](const string& target) {
        unordered_map<string, vector<int>> index;
        if (target.length() >= k) {
            for (int i = 0; i <= static_cast<int>(target.length()) - k; ++i) {
                string kmer = target.substr(i, k);
                index[kmer].push_back(i);
            }
        }
        return index;
    };

    unordered_map<string, vector<int>> kmer_index_seq = build_kmer_index(seq);
    unordered_map<string, vector<int>> kmer_index_rc = build_kmer_index(rc_seq);

    // 3. 锚点比对函数（k-mer优化版）
    auto align_anchor = [&](const string& anchor, const string& target, 
                            const unordered_map<string, vector<int>>& index) -> pair<int, int> {
        const int L = anchor.length();
        const int max_ed = error_count;
        unordered_set<int> candidates;
        const int step = max_ed + 1;  // 鸽巢原理步长

        // 3.1 收集候选位置
        for (int i = 0; i <= L - k; i += step) {
            string kmer = anchor.substr(i, k);
            if (index.find(kmer) != index.end()) {
                for (int pos : index.at(kmer)) {
                    int start = pos - i;
                    if (start >= 0 && start <= static_cast<int>(target.length()) - L) {
                        candidates.insert(start);
                    }
                }
            }
        }

        // 3.2 候选位置比对
        int min_dist = INT_MAX;
        int best_end = -1;
        for (int pos : candidates) {
            string window = target.substr(pos, L);
            int dist = edlib_align(anchor, window, max_ed);
            if (dist != -1 && dist <= max_ed && dist < min_dist) {
                min_dist = dist;
                best_end = pos + L - 1;  // 记录结束位置
            }
        }
        return {min_dist, best_end};
    };

    // 4. 执行四组比对（正向/反向 x 左锚点/右锚点）
    auto [left_dist, left_end] = align_anchor(anchors[0], seq, kmer_index_seq);
    auto [right_dist, right_start] = align_anchor(anchors[1], seq, kmer_index_seq);
    auto [left_rc_dist, left_rc_end] = align_anchor(anchors[0], rc_seq, kmer_index_rc);
    auto [right_rc_dist, right_rc_start] = align_anchor(anchors[1], rc_seq, kmer_index_rc);

    // 5. 判断命中逻辑（同原函数）
    int minscore = anchors[0].size() - error_count;
    bool fwd_hit = (left_dist <= error_count || right_dist <= error_count) && 
                  (left_rc_dist > error_count && right_rc_dist > error_count);
    bool rev_hit = (left_rc_dist <= error_count || right_rc_dist <= error_count) && 
                  (left_dist > error_count && right_dist > error_count);

    // 6. 提取CID
    if (fwd_hit || rev_hit) {
        const string& tureseq = fwd_hit ? seq : rc_seq;
        int leftPosE = fwd_hit ? left_end : left_rc_end;
        int rightPosS = fwd_hit ? right_start : right_rc_start;

        if (left_dist <= error_count && right_dist <= error_count) {
            int width = rightPosS - leftPosE;
            if (width >= widths[0] && width <= widths[1]) {
                cid = tureseq.substr(leftPosE, width);
                align_status = (fwd_hit ? "+\t" : "-\t") + to_string(leftPosE) + "\t2_anchor";
            }
        } else if (right_dist <= error_count) {
            cid = tureseq.substr(rightPosS - widths[2], widths[2]);
            align_status = (fwd_hit ? "+\t" : "-\t") + to_string(rightPosS - widths[2]) + "\tright_anchor";
        } else if (left_dist <= error_count) {
            cid = tureseq.substr(leftPosE, widths[2]);
            align_status = (fwd_hit ? "+\t" : "-\t") + to_string(leftPosE) + "\tleft_anchor";
        }
    }

    return FsUtils::reverse_complement(cid);
}

std::string CidExtractPairwiseAlign(std::string sequence, Config& config, string& align_status) {
    string cid = "";
    const auto& anchors = config.anchor_seqs;
    const auto& widths = config.width_params;
    const int error_count = config.error_threshold;

    const int sub_len = std::min(150, static_cast<int>(sequence.size()));
    const std::string seq = sequence.substr(sequence.size()-sub_len, sub_len);
    const std::string rc_seq = FsUtils::reverse_complement(sequence).substr(sequence.size()-sub_len, sub_len);

    auto left_dna5 = anchors[0] | seqan3::views::char_to<seqan3::dna5>;
    auto right_dna5 = anchors[1] | seqan3::views::char_to<seqan3::dna5>;
    auto method = seqan3::align_cfg::method_local{};
    seqan3::align_cfg::scoring_scheme scheme{ seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}} };
    seqan3::align_cfg::gap_cost_affine gap_costs{ seqan3::align_cfg::open_score{0}, seqan3::align_cfg::extension_score{-1} };
    auto aligncfg = method | scheme | gap_costs;    
    auto seq_dna5 = seq | seqan3::views::char_to<seqan3::dna5>;
    auto seq_rc_dna5 = rc_seq | seqan3::views::char_to<seqan3::dna5>;

    auto left_results = seqan3::align_pairwise(std::tie(left_dna5, seq_dna5), aligncfg);
    auto right_results = seqan3::align_pairwise(std::tie(right_dna5, seq_dna5), aligncfg);
    auto left_rc_results = seqan3::align_pairwise(std::tie(left_dna5, seq_rc_dna5), aligncfg);
    auto right_rc_results = seqan3::align_pairwise(std::tie(right_dna5, seq_rc_dna5), aligncfg);
    auto& left_res = *left_results.begin(); auto& right_res = *right_results.begin(); auto& left_rc_res = *left_rc_results.begin(); auto& right_rc_res = *right_rc_results.begin();
    uint32_t leftPosEP = left_res.sequence2_end_position(); uint32_t rightPosSP = right_res.sequence2_begin_position();
    uint32_t leftPosEM = left_rc_res.sequence2_end_position(); uint32_t rightPosSM = right_rc_res.sequence2_begin_position();
    uint32_t leftPosE, rightPosS;
    string tureseq;
    int width, hit = 0;

    int minscore1 = anchors[0].size() - error_count;
    int minscore2 = anchors[1].size() - error_count;
    //cout << "left score: " << left_res.score() << "right score: " << right_res.score() << "minscore " << minscore << endl;
    //cout << leftPosEP << "\t" << rightPosSP << "\t" << leftPosEM << "\t" << rightPosSM << endl;
    int leftscore, rightscore;
    if ((left_res.score() >= minscore1 || right_res.score() >= minscore2) && (left_rc_res.score() < minscore1 && right_rc_res.score() < minscore2)) {
        hit++;
        leftPosE = leftPosEP;
        rightPosS = rightPosSP;
        tureseq = seq;
        leftscore = left_res.score();
        rightscore = right_res.score();
        align_status = "+\t";
    }
    else if ((left_res.score() < minscore1 && right_res.score() < minscore2) && (left_rc_res.score() >= minscore1 || right_rc_res.score() >= minscore2)) {
        hit++;
        leftPosE = leftPosEM;
        rightPosS = rightPosSM;
        tureseq = rc_seq;
        leftscore = left_rc_res.score();
        rightscore = right_rc_res.score();
        align_status = "-\t";
    }
    if (hit > 0) {
        if (leftscore >= minscore1 && rightscore >= minscore2) {
            width = rightPosS - leftPosE;
            if (width >= widths[0] && width <= widths[1]) {
                cid = tureseq.substr(leftPosE, width);
                align_status = align_status + to_string(leftPosE) + "\t" + "2_anchor";
            }
        }
        else if (rightscore >= minscore2 && rightPosS >= widths[2]) {
            cid = tureseq.substr(rightPosS - widths[2], widths[2]);
            align_status = align_status + to_string(rightPosS - widths[2]) + "\t" + "right_anchor";
        }
        else if (leftscore >= minscore1) {
            cid = tureseq.substr(leftPosE, widths[2]);
            align_status = align_status + to_string(leftPosE) + "\t" + "left_anchor";
        }
    }

    cid = FsUtils::reverse_complement(cid); // CID whitelist is 3->5
    return(cid);
}
vector<streampos> get_fastq_offsets(const string& fastq_path) {
    ifstream inFile(fastq_path);
    if (!inFile) {
        cerr << "Error opening file: " << fastq_path << endl;
        return {};
    }

    vector<streampos> offsets;
    streampos pos = inFile.tellg();
    string line;

    while (true) {
        // 阶段1：寻找记录起始标识
        while (getline(inFile, line)) {
            if (!line.empty() && line[0] == '@') {
                offsets.push_back(pos);
                break;
            }
            pos = inFile.tellg();
        }
        if (inFile.eof()) break;

        // 阶段2：跳过完整记录
        string seq_line, plus_line, qual_line;
        if (!getline(inFile, seq_line) || 
            !getline(inFile, plus_line) || 
            !getline(inFile, qual_line)) {
            break;
        }
        
        // 更新位置
        pos = inFile.tellg();
    }

    return offsets;
}

struct FastqRecord {
    string id;
    string seq;
};

void handleFastqChunk(const string& fastqFile, 
                     const vector<streampos>& read_offsets,
                     uint64_t start_idx, uint64_t end_idx,
                     Config& config, const string& tmpOutputPath, 
                     uint64_t& local_processed, uint64_t& local_cid_count) {
    
    ifstream inFile(fastqFile);
    if (!inFile) {
        cerr << "Error opening file: " << fastqFile << endl;
        return;
    }

    ofstream outFile(tmpOutputPath);
    if (!outFile) {
        cerr << "Error opening temp output: " << tmpOutputPath << endl;
        return;
    }

    // 优化IO性能
    ios_base::sync_with_stdio(false);
    inFile.tie(nullptr);

    local_processed = 0;
    local_cid_count = 0;

    // 处理当前分块内的所有记录
    for (uint64_t idx = start_idx; idx < end_idx; ++idx) {
        inFile.clear(); // 清除可能的EOF标志
        inFile.seekg(read_offsets[idx]);
        
        string header, sequence, plus, quality;
        if (!getline(inFile, header)) break;
        if (!getline(inFile, sequence)) break;
        if (!getline(inFile, plus)) break;
        if (!getline(inFile, quality)) break;

        // 处理CID提取
        string alignstat;
        string cid = CidExtractPairwiseAlign(sequence, config, alignstat);

        // 写入临时结果文件
        if (!cid.empty()) {
            // 移除@符号作为ID
            string record_id = header.substr(1);
            outFile << record_id << '\t'
                    << alignstat << '\t'
                    << cid << '\n';
            ++local_cid_count;
        }

        ++local_processed;

        // // 定期输出进度
        // if (local_processed % 10000 == 0) {
        //     cout << "Thread [" << this_thread::get_id() 
        //          << "] processed " << local_processed 
        //          << " records" << endl;
        // }
    }

    // 显式关闭文件
    inFile.close();
    outFile.close();
}

// 修改：合并临时文件到最终输出
void merge_temp_files(const vector<string>& tmp_paths, const string& output_path) {
    ofstream final_out(output_path);
    for (const auto& tmp_path : tmp_paths) {
        ifstream in_tmp(tmp_path);
        if (in_tmp) {
            final_out << in_tmp.rdbuf();
            in_tmp.close();
            remove(tmp_path.c_str()); // 删除临时文件
        }
    }
}

string extractcid_fastq(const string& fastq_path,
                   const string& output_path,
                   Config& config) {
    // ===== 1. 预计算文件偏移量 =====
    cout << "Precomputing FASTQ record offsets..." << endl;
    vector<streampos> read_offsets = get_fastq_offsets(fastq_path);
    if (read_offsets.empty()) {
        return "Error: Failed to compute record offsets";
    }
    
    size_t total_records = read_offsets.size();
    cout << "Found " << total_records << " records in " << fastq_path << endl;

    // ===== 2. 任务分块 =====
    vector<pair<size_t, size_t>> chunks;
    const size_t chunk_size = total_records / config.num_threads;
    size_t remaining = total_records % config.num_threads;
    size_t start_idx = 0;
    
    for (size_t i = 0; i < config.num_threads; ++i) {
        size_t end_idx = start_idx + chunk_size + (i < remaining ? 1 : 0);
        chunks.emplace_back(start_idx, end_idx);
        start_idx = end_idx;
    }

    // ===== 3. 多线程处理 =====
    vector<thread> workers;
    vector<string> tmp_outputs;
    vector<uint64_t> chunk_processed(config.num_threads, 0);
    vector<uint64_t> chunk_cid_count(config.num_threads, 0);

    // 创建临时文件路径
    for (size_t tid = 0; tid < config.num_threads; ++tid) {
        tmp_outputs.push_back(output_path + ".part" + to_string(tid) + ".tmp");
    }

    cout << "Starting parallel processing with " << config.num_threads << " threads..." << endl;
    for (size_t tid = 0; tid < config.num_threads; ++tid) {
        workers.emplace_back(
            handleFastqChunk, 
            cref(fastq_path),
            cref(read_offsets),
            chunks[tid].first,
            chunks[tid].second,
            ref(config),
            cref(tmp_outputs[tid]),
            ref(chunk_processed[tid]),
            ref(chunk_cid_count[tid])
        );
    }

    // ===== 4. 等待线程完成 =====
    for (auto& t : workers) {
        t.join();
    }

    // ===== 5. 合并结果 =====
    cout << "Merging temporary files..." << endl;
    merge_temp_files(tmp_outputs, output_path);

    // ===== 6. 汇总统计 =====
    total_processed = accumulate(chunk_processed.begin(), chunk_processed.end(), 0ULL);
    valid_cid_count = accumulate(chunk_cid_count.begin(), chunk_cid_count.end(), 0ULL);

    string summary = "Total processed sequences: " + to_string(total_processed) 
                   + "\tValid CID counts: " + to_string(valid_cid_count) + "\n";
    
    // 释放偏移量向量内存
    vector<streampos>().swap(read_offsets);
    
    return summary;
}

} // namespace CIDExtract