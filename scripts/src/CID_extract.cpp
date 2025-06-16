// cidmap_extract.cpp
#include "CID_extract.h"
#include "utils.h"
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>
#include <memory>
#include <iostream>
#include <string>
#include <unordered_map>
#include <ranges>
#include <htslib/sam.h>
#include <htslib/hts.h>

using namespace std;

namespace CIDExtract {

using namespace seqan3::literals;

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
void handleBamChunk(const string& bamFile, uint64_t start, uint64_t end, GeneRanges & geneRanges, 
    Config& config,
     ofstream& outputFile, mutex& mtx,
     uint64_t& total_processed,
     uint64_t& valid_cid_count,
     mutex& count_mtx) {
    samFile* fp = sam_open(bamFile.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(fp);
    bam1_t* record = bam_init1();
    
    uint64_t local_processed = 0;
    uint64_t local_cid_count = 0;

    uint64_t count = 0;
    while (sam_read1(fp, header, record) >= 0) {
        if (record->core.flag & BAM_FUNMAP) {
            continue;
        }

        if (++count < start) {
            continue;
        }

        if (count > end) {
            break;
        }
        if (count % 10000 == 0) {
            cout << count << endl;
        }
        local_processed++;
        ReadInfo read;
        read.id = bam_get_qname(record);

        uint8_t* seq = bam_get_seq(record);
        std::string seq_str;

        for (int i = 0; i < record->core.l_qseq; ++i) {
            seq_str.push_back(seq_nt16_str[bam_seqi(seq, i)]);
        }
        read.seq = seq_str;

        int tid = record->core.tid;
        const char* chrName = header->target_name[tid];
        string chrNameStr(chrName);
        uint64_t pos = record->core.pos;
        uint64_t end = bam_endpos(record);
        //cout << chrNameStr << "\t" << pos << "\t" << end << endl;
        read.gene = "none";
        for (const auto& gene : geneRanges[chrNameStr]) {
            if (gene.first.empty()) {
                continue;
            }
            if (pos <= gene.second.second && end >= gene.second.first) {
                //cout << gene.first << "\t" << gene.second.first << "\t" << gene.second.second << endl;
                if (read.gene == "none") {
                    read.gene = gene.first;
                }
                else {
                    read.gene = read.gene + "," + gene.first;
                }

            }
        }

        //extract cid
        string alignstat;
        string cid = CidExtractPairwiseAlign(read.seq, config, alignstat);
        //cout << cid << endl;
        
        if (cid != "") {
            local_cid_count++;
            mtx.lock();
            // Lock the output file stream to prevent multiple threads writing to it simultaneously
            outputFile << read.id << '\t' << alignstat << "\t" << cid << "\t" << read.gene << "\t" << pos << "\t" << end << endl;
        }
        mtx.unlock();
    }

    sam_close(fp);
    bam_destroy1(record);
    bam_hdr_destroy(header);

    std::lock_guard<std::mutex> lock(count_mtx);
    total_processed += local_processed;
    valid_cid_count += local_cid_count;
}


// 辅助函数：分割字符串
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        if (!token.empty()) {
            tokens.push_back(token);
        }
    }
    return tokens;
}

void load_gene_ranges(const std::string& gtf_path, GeneRanges& gene_ranges) {
    std::ifstream file(gtf_path);
    std::string line;
    
    while (std::getline(file, line)) {
        // 跳过注释行
        if (line[0] == '#') continue;

        std::vector<std::string> fields = split(line, '\t');
        if (fields.size() < 9) continue;

        // 解析基础字段
        const std::string& seqname = fields[0];
        const std::string& feature_type = fields[2];
        uint64_t start = std::stoull(fields[3]);
        uint64_t end = std::stoull(fields[4]);

        // 只处理gene类型记录[6,9,11](@ref)
        if (feature_type != "gene") continue;

        // 解析属性字段
        std::string gene_id, gene_name;
        std::vector<std::string> attributes = split(fields[8], ';');
        for (const auto& attr : attributes) {
            auto key_value = split(attr, ' ');
            if (key_value.size() < 2) continue;
            
            if (key_value[0] == "gene_id") {
                gene_id = key_value[1].substr(1, key_value[1].size()-2); // 去除引号
            } else if (key_value[0] == "gene_name") {
                gene_name = key_value[1].substr(1, key_value[1].size()-2);
            }
        }

        // 优先使用gene_name，若无则使用gene_id[6,9](@ref)
        std::string gene_identifier = gene_name.empty() ? gene_id : gene_name;
        if (!gene_identifier.empty()) {
            // 更新基因范围（考虑同基因多记录的情况）[6](@ref)
            auto& gene_map = gene_ranges[seqname];
            if (gene_map.find(gene_identifier) == gene_map.end()) {
                gene_map[gene_identifier] = {start, end};
            } else {
                // 合并重叠区域，取最小start和最大end[6](@ref)
                auto& current = gene_map[gene_identifier];
                current.first = std::min(current.first, start);
                current.second = std::max(current.second, end);
            }
        }
    }
}

uint64_t count_bam_records(const std::string& bam_path) {
    htsFile* fp = sam_open(bam_path.c_str(), "r");
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> 
        header(sam_hdr_read(fp), bam_hdr_destroy);
    bam1_t* rec = bam_init1();
    uint64_t count = 0;
    
    while (sam_read1(fp, header.get(), rec) >= 0) ++count;
    
    bam_destroy1(rec);
    sam_close(fp);
    return count;
}

string process_bam(const std::string& bam_path,
                const std::string& output_path,
                GeneRanges& gene_ranges,
                Config& config) {
    std::ofstream output(output_path);
    std::mutex mtx;
    uint64_t numThreads = config.num_threads;
    const uint64_t total_records = count_bam_records(bam_path);
    const uint64_t chunk_size = total_records / numThreads;
    uint64_t remainder = total_records % numThreads;
    
    vector<thread> threads;
    std::mutex count_mtx;  // 新增互斥锁
    uint64_t total_processed = 0;
    uint64_t valid_cid_count = 0;

    uint64_t start = 0, end = 0;
    for (uint64_t i = 0; i < numThreads; i++) {
        start = end;
        end = start + chunk_size;
        if (i == numThreads - 1) {
            end += remainder;
        }
        cout << start << "-" << end << endl;
        threads.emplace_back(handleBamChunk, bam_path, start, end, ref(gene_ranges), ref(config), ref(output), ref(mtx),
        ref(total_processed),ref(valid_cid_count),ref(count_mtx));
    }

    for (auto& t : threads) {
        t.join();
    }
    string summary = "Total mapped reads(segments) within gene body: " + to_string(total_processed) + "\t" + "Reads with valid CID: " + to_string(valid_cid_count) + "\n";
    return summary;
}

} // namespace CIDExtract