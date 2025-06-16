// CIDindex.cpp
#include "CID_index.h"
#include <fstream>
#include <thread>
#include <sstream>
#include <vector>
#include <mutex>
#include <bitset>
#include <set>
#include <map>
#include <cereal/archives/binary.hpp>
using namespace std;

namespace CIDindex {

size_t KMERLEN = 6;
size_t BUCKETNUM = 6;
const int K = 25; 
const int BITS_PER_BASE = 2; 

struct Seqinfo {
    bitset<K* BITS_PER_BASE>  encode_seqs;
    uint32_t  xpos;
    uint32_t  ypos;
    template<class Archive>
    void serialize(Archive& ar) {
        ar(encode_seqs, xpos, ypos);
    }
};

std::size_t getCurrentMemoryUsage() {
    std::ifstream statm("/proc/self/statm");
    if (!statm) {
        std::cerr << "Failed to open /proc/self/statm" << std::endl;
        return 0;
    }

    std::string line;
    std::getline(statm, line);
    std::istringstream iss(line);
    std::size_t pages;
    iss >> pages;

    const std::size_t page_size = sysconf(_SC_PAGESIZE);
    return pages * page_size;
}

uint32_t murmurHash3(string seq, uint32_t seed) {
    uint32_t out;
    MurmurHash3_x64_128(seq.data(), seq.size(), seed, &out);
    return out;
}

uint32_t min_hash(string seq, size_t len, uint32_t seed, uint8_t &murpos) {
    uint64_t min_val = UINT64_MAX;
    size_t min_index = 0;

    for (uint8_t i = 0; i <= (seq.size() - len); i++) {
        string seqsub = seq.substr(i, len);
        uint32_t hasht = murmurHash3(seqsub, seed);

        if (hasht < min_val) {
            min_val = hasht;
            min_index = i;
        }
    }

    murpos = static_cast<uint8_t>(min_index);
    return min_val;
}

void bucket_build(map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>& bucketi, string seq, uint32_t hashTotal) {
    uint8_t murpos;
    for (uint8_t s = 0; s < BUCKETNUM; s++) {
        uint32_t hashMin = min_hash(seq, KMERLEN, s, murpos);
        bucketi[s][hashMin].push_back(make_pair(hashTotal, murpos));
    }
}

bitset<K* BITS_PER_BASE> encode_sequence(const string& sequence) {
    bitset<K* BITS_PER_BASE> encoded_seq;
    int bit_index = 0;
    for (char base : sequence) {
        switch (base) {
        case 'A': encoded_seq.set(bit_index, 0); encoded_seq.set(bit_index + 1, 0); break;
        case 'T': encoded_seq.set(bit_index, 0); encoded_seq.set(bit_index + 1, 1); break;
        case 'C': encoded_seq.set(bit_index, 1); encoded_seq.set(bit_index + 1, 0); break;
        case 'G': encoded_seq.set(bit_index, 1); encoded_seq.set(bit_index + 1, 1); break;
        }
        bit_index += BITS_PER_BASE;
    }
    return encoded_seq;
}

void gene_index_precess_chunk(string tableFile, uint64_t start_line, uint64_t end_line, 
    map<string, map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>>& TotalBuckets,
    map <string, map<uint32_t, Seqinfo>>& seqBuckets,
    mutex& mtx) {
    
    ifstream inFile(tableFile);
    string line;
    uint64_t current_line = 0;
    set<uint32_t> uniqueSet;

    while (getline(inFile, line)) {
        if (current_line < start_line) {
            current_line++;
            continue;
        }
        if (current_line > end_line) {
            break;
        }
        current_line++;

        vector<string> columns;
        stringstream ss(line);
        string column;
        while (getline(ss, column, '\t')) {
            columns.push_back(column);
        }

        if (columns.size() < 4) {
            cerr << "Invalid line format: " << line << endl;
            continue;
        }

        string cidseq = columns[0];
        uint32_t xpos, ypos;
        try {
            xpos = stoul(columns[1]);
            ypos = stoul(columns[2]);
        } catch (const exception& e) {
            cerr << "Error parsing coordinates: " << line << endl;
            continue;
        }

        string gene_str = columns[3];
        vector<string> genes;
        stringstream gene_ss(gene_str);
        string gene;
        while (getline(gene_ss, gene, ';')) {
            genes.push_back(gene);
        }

        uint32_t hashTotal = murmurHash3(cidseq, 0);
        if (uniqueSet.count(hashTotal) == 0) {
            uniqueSet.insert(hashTotal);

            map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>> bucketi;
            bucket_build(bucketi, cidseq, hashTotal);

            Seqinfo seqinfoi;
            seqinfoi.encode_seqs = encode_sequence(cidseq);
            seqinfoi.xpos = xpos;
            seqinfoi.ypos = ypos;

            mtx.lock();
            for (const string& gene : genes) {
                for (auto& [s, inner_map] : bucketi) {
                    for (auto& [hashMin, vec] : inner_map) {
                        TotalBuckets[gene][s][hashMin].insert(
                            TotalBuckets[gene][s][hashMin].end(), 
                            vec.begin(), vec.end()
                        );
                    }
                }
                seqBuckets[gene][hashTotal] = seqinfoi;
            }
            mtx.unlock();
        }
    }
}

void gene_index_process(string tableFile, int numThreads,
    map<string, map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>>& TotalBuckets,
    map <string, map<uint32_t, Seqinfo>>& seqBuckets) {
    
    mutex mtx;
    ifstream inFile(tableFile);
    uint64_t total_lines = count(istreambuf_iterator<char>(inFile), 
                                istreambuf_iterator<char>(), '\n');
    inFile.close();

    total_lines++;  // 修正行数统计
    uint64_t chunk_size = total_lines / numThreads;
    vector<thread> threads;

    for (int i = 0; i < numThreads; ++i) {
        uint64_t start = i * chunk_size;
        uint64_t end = (i == numThreads-1) ? total_lines : start + chunk_size;
        threads.emplace_back(gene_index_precess_chunk, 
                            tableFile, start, end,
                            ref(TotalBuckets), ref(seqBuckets), ref(mtx));
    }

    for (auto& t : threads) {
        t.join();
    }
}
int calculate_edit_threshold(int N) {
    const int max_d = 10; // 假设最大可能d不超过10
    for (int d = 1; d <= max_d; ++d) {
        double p_error = 0.0;
        for (int k = 0; k <= d; ++k) {
            // 计算二项分布概率项
            double comb = std::tgamma(25 + 1) / (std::tgamma(k + 1) * std::tgamma(25 - k + 1));
            p_error += comb * std::pow(0.75, k) * std::pow(0.25, 25 - k);
        }
        if (N * p_error > 0.001) {
            return d - 1;
        }
    }
    return max_d; // 若未找到，返回安全值
}

void build_index(string tableFile, int kmerlen, int bucketnum, int numThreads, string outFold) {
    KMERLEN = kmerlen;
    BUCKETNUM = bucketnum;
    
    map<string, map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>> TotalBuckets;
    map <string, map<uint32_t, Seqinfo>> seqBuckets;

    cout << "Building index from table..." << endl;
    gene_index_process(tableFile, numThreads, TotalBuckets, seqBuckets);

    string outFm = outFold + ".bin";
    ofstream os(outFm, ios::binary);
    cereal::BinaryOutputArchive oarchive(os);
    oarchive(TotalBuckets, seqBuckets);
    os.close();

    cout << "Index built successfully. Output saved to " << outFm << endl;

    // 生成阈值表格文件
    string thresholdFile = outFold + "_thresholds.tsv";
    ofstream thresholdStream(thresholdFile);
    if (!thresholdStream) {
        cerr << "Error opening threshold table file: " << thresholdFile << endl;
        return;
    }
    // 写入表头
    thresholdStream << "Gene\tCID_Count\tThreshold\n";
    // 遍历每个基因计算阈值
    for (const auto& [gene, cidMap] : seqBuckets) {
        size_t cid_count = cidMap.size();
        // 这里替换为实际的阈值计算公式
        int threshold = 0;
        if (cid_count > 0) {
            threshold = calculate_edit_threshold(cid_count);
        }
        thresholdStream << gene << "\t" << cid_count << "\t" << threshold << "\n";
    }

    thresholdStream.close();
    cout << "Threshold table saved to " << thresholdFile << endl;
}

}