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
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/unordered_map.hpp>
using namespace std;

namespace CIDindexPrecise {
constexpr size_t SHARD_COUNT = 64;
size_t KMERLEN = 6;
size_t BUCKETNUM = 6;
const int K = 25; 
const int BITS_PER_BASE = 2; 

#pragma pack(push, 1)
struct SeqInfo {
    uint32_t xpos;        // 使用位域压缩坐标（支持131071坐标范围）
    uint32_t ypos;
    template<class Archive>
    void serialize(Archive& ar) {
        ar(xpos, ypos);
    }
};
#pragma pack(pop)

struct Bucket {
    unordered_map<uint64_t, SeqInfo> sub_buckets; // murpos作为二级键
    template<class Archive>
    void serialize(Archive& ar) {
        ar(sub_buckets); // 新增序列化函数
    }
};
std::array<std::mutex, SHARD_COUNT> shard_mutexes{};
array<unordered_map<uint32_t, Bucket>, SHARD_COUNT> sharded_table;

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
    MurmurHash3_x86_32(seq.data(), seq.size(), seed, &out);
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


uint64_t encode_sequence(const string& sequence) {
    uint64_t encoded = 0;
    for (char base : sequence) {
        encoded <<= 2;
        switch (toupper(base)) { // 支持小写字母转换
            case 'A': encoded |= 0b00; break;
            case 'T': encoded |= 0b01; break;
            case 'C': encoded |= 0b10; break;
            case 'G': encoded |= 0b11; break;
            default: throw invalid_argument("Invalid base: " + string(1, base));
        }
    }
    return encoded;
}

void index_process_chunk(string tableFile, uint64_t start_line, uint64_t end_line) {
    
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

        // if (columns.size() < 4) {
        //     cerr << "Invalid line format: " << line << endl;
        //     continue;
        // }

        string cidseq = columns[0];
        uint32_t xpos, ypos;
        try {
            xpos = stoul(columns[1]);
            ypos = stoul(columns[2]);
        } catch (const exception& e) {
            cerr << "Error parsing coordinates: " << line << endl;
            continue;
        }

        //uint64_t encode_seqs = encode_sequence(cidseq);
        //uint32_t hashTotal = murmurHash3(cidseq, 0);

        SeqInfo seqinfoi;
        auto encoded_seq = encode_sequence(cidseq);
        seqinfoi.xpos = xpos;
        seqinfoi.ypos = ypos;
        uint8_t murpos;
        uint32_t hashMin = min_hash(cidseq, KMERLEN, 0, murpos);

        // 分片存储
        size_t shard_idx = (hashMin * 11400714819323198485ULL) >> 58; // Fibonacci哈希
        lock_guard<mutex> lock(shard_mutexes[shard_idx]);
        sharded_table[shard_idx][hashMin].sub_buckets[encoded_seq] = seqinfoi;
    }

}

void build_index_precise(string tableFile, int kmerlen, int numThreads, string outFold) {
    KMERLEN = kmerlen;

    cout << "Building index from table..." << endl;
    ifstream inFile(tableFile);
    uint64_t total_lines = count(istreambuf_iterator<char>(inFile), 
                                istreambuf_iterator<char>(), '\n');
    inFile.close();
    vector<thread> workers;
    uint64_t chunk_size = total_lines / numThreads;
    
    for (int i = 0; i < numThreads; ++i) {
        workers.emplace_back([=] {
            index_process_chunk(tableFile, i*chunk_size, (i+1)*chunk_size);
        });
    }
    for (auto& t : workers) t.join();

    string outFm = outFold + ".precise.bin";

    ofstream os(outFm, ios::binary);
    cereal::BinaryOutputArchive oarchive(os);
    oarchive(sharded_table);
    os.close();

    cout << "Index built successfully. Output saved to " << outFm << endl;
}

}