#include "CID_mapping_precise.h"
#include <chrono>
#include <thread>
#include <unordered_set>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/unordered_map.hpp>

using namespace std;
namespace CIDmappingPrecise {
constexpr size_t SHARD_COUNT = 64;
size_t KMERLEN = 6;
size_t BUCKETNUM = 6;
const int K = 25;
const int BITS_PER_BASE = 2;
const int BASES_PER_BYTE = 4;
size_t hitcount = 0;

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
array<unordered_map<uint32_t, Bucket>, SHARD_COUNT> sharded_table;

struct ReadInfo {
    string id;
    string seq;
    pair<uint64_t, uint64_t> pos;
    string gene;
    string strand;
};

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


void kmer_hash(string seq, size_t len, uint64_t seed, vector<uint64_t>& hashvec) {
    for (size_t i = 0; i <= (seq.size() - len); i++)
    {
        string seqsub = seq.substr(i, len);
        uint64_t hasht = murmurHash3(seqsub, seed);
        hashvec.push_back(hasht);
    }
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
string decode_sequence(uint64_t encoded_seq) {
    string decoded;
    decoded.reserve(25);  // 预分配25个字符空间
    
    // 从最高有效位开始解码
    uint64_t mask = 0xC000000000000000;  // 二进制11000000...（用于提取最高2位）
    
    for (int i = 0; i < 25; ++i) {
        // 提取当前最高2位
        uint8_t bits = (encoded_seq & mask) >> (62 - i*2);
        
        // 转换为碱基字符
        switch(bits) {
            case 0b00: decoded += 'A'; break;
            case 0b01: decoded += 'T'; break;
            case 0b10: decoded += 'C'; break;
            case 0b11: decoded += 'G'; break;
        }
        
        // 更新掩码处理下个2位
        mask >>= 2;
    }
    
    return decoded;
}

vector<string> generate_candidates(const string& query) {
    vector<string> candidates;
    candidates.reserve(129);

    // 替换突变
    for (int i = 0; i < query.size(); ++i) {
        string temp = query;
        for (char c : {'A','T','C','G'}) {
            if (c != temp[i]) {
                temp[i] = c;
                candidates.emplace_back(temp);
            }
        }
    }

    // 插入/删除突变
    if(query.size() == 24) {
        for (int i = 0; i <= query.size(); ++i) {
            for (char c : {'A','T','C','G'}) {
                string temp = query;
                temp.insert(i, 1, c);
                candidates.emplace_back(temp);
            }
        }
    } else if(query.size() == 26) {
        for (int i = 0; i < query.size(); ++i) {
            string temp = query;
            temp.erase(i, 1);
            candidates.emplace_back(temp);
        }
    }

    return candidates;
}

void loadQueryReads(string readFile, vector<ReadInfo>& readsInfols) {
    ifstream inFile(readFile);
    if (!inFile) {
        cerr << "Error: failed to open " << readFile << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    while (getline(inFile, line)) {
        istringstream ss(line);
        string readid, strand, pos, status, seq;
        getline(ss, readid, '\t');
        getline(ss, strand, '\t');
        getline(ss, pos, '\t');
        getline(ss, status, '\t');
        getline(ss, seq, '\t');
        ReadInfo readInfo;
        readInfo.id = readid;
        readInfo.seq = seq;
        readInfo.strand = "*";
        readsInfols.push_back(readInfo);
    }
}


void matchCidChunk(vector<ReadInfo>& readsInfols,
    ofstream& outputFile, int start, int end, mutex& mtx) {
    vector<string> outputBuffer;
    const size_t bufferSize = 1000; 
    size_t hitcounti = 0;
    for (int i = start; i <= end; i++) {
        string seqi = readsInfols[i].seq;
        string strand = readsInfols[i].strand;
        string readid = readsInfols[i].id + "\t" + readsInfols[i].seq;
        string resPos, resSeq;

        if(seqi.size() < 20 || seqi.size() > 30){
            continue;
        }
        vector<pair<uint32_t, uint32_t>> results;
        uint8_t murpos;
        uint32_t hashMin = min_hash(seqi, KMERLEN, 0, murpos);
        //cout << seqi << " minHash :" << hashMin << endl;
        uint64_t seqi_ed = encode_sequence(seqi);
        size_t shard_idx = (hashMin * 11400714819323198485ULL) >> 58;
        uint64_t candidate_ed;
        mtx.lock();
        auto& shard = sharded_table[shard_idx];
        mtx.unlock();
        //cout << "shard size: " << shard.size() << hashMin << endl;
        auto it = shard.find(hashMin);
        if(it != shard.end()){
            //cout << seqi << " found minhash" << endl;
            auto& sub_buckets = it->second.sub_buckets;
            auto sub_it = sub_buckets.find(seqi_ed);
            if(sub_it != sub_buckets.end()){
                auto &seqInfo = sub_it->second;
                resPos = to_string(seqInfo.xpos) + "_" + to_string(seqInfo.ypos);
                outputBuffer.emplace_back(readid + "\t" + resPos + "\t" + seqi + "\t0\t" + "*" + "\t" + 
                "*" + "\t" + strand);
                hitcounti++;
                continue;
            }
        }

        auto candidates = generate_candidates(seqi);
        for (const auto& candidate : candidates) {
            hashMin = min_hash(candidate, KMERLEN, 0, murpos);
            candidate_ed = encode_sequence(candidate);
            shard_idx = (hashMin * 11400714819323198485ULL) >> 58;
            if (shard_idx >= SHARD_COUNT) continue;
            mtx.lock();
            auto& shard = sharded_table[shard_idx];
            mtx.unlock();
            auto it = shard.find(hashMin);
            if (it == shard.end()) continue;
            //cout << candidate << " found minhash with error" << endl;
            auto& sub_buckets = it->second.sub_buckets;
            auto sub_it = sub_buckets.find(candidate_ed);
            if (sub_it == sub_buckets.end()){
                continue;
            }else{
                auto &seqInfo = sub_it->second;
                resPos = to_string(seqInfo.xpos) + "_" + to_string(seqInfo.ypos);
                outputBuffer.emplace_back(readid + "\t" + resPos + "\t" + candidate + "\t1\t" + "*" + "\t" + 
                "*" + "\t" + strand);
                hitcounti++;
                break;
            }
        }
    }

    mtx.lock();
    hitcount = hitcount+hitcounti;
    for (const auto& line : outputBuffer) {
        outputFile << line << endl;
    }
    mtx.unlock();

}

string matchCidThread(vector<ReadInfo>& readsInfols,
    string outputFileName, int numThreads) {

    vector<thread> threads;
    mutex mtx;
    int n_records = readsInfols.size();
    int chunk_size = n_records / numThreads;
    int remainder = n_records % numThreads;

    ofstream outputFile(outputFileName);
    outputFile << "readid\tquerySeq\tcidPos\trefSeq\teditDi\teditDi_ref\tgene\tstrand" << endl;
    int start = 0, end = 0;
    for (int i = 0; i < numThreads; i++) {
        start = end;
        end = start + chunk_size;
        if (i == numThreads - 1) {
            end += remainder - 1;
        }
        cout << start << "-" << end << endl;
        threads.emplace_back(matchCidChunk, ref(readsInfols),ref(outputFile), start, end, ref(mtx));
    }

    for (auto& t : threads) {
        t.join();
    }
    
    double ratio = hitcount/n_records * 100;
    string summary = "Total sequences: " + to_string(n_records) + "\t" + "Passed sequences: " + to_string(hitcount) + "\tPassed ratio: " + to_string(ratio) + "%\n";
    return summary;
}

string CIDmap_processer_precise(string readFile, string indexFold, int thread, int kmerlen, string outFile) {
    KMERLEN = kmerlen;
    vector<ReadInfo> readsInfols;

    // load dp index
    cout << "loading index" << endl;
    ifstream is(indexFold, ios::binary);
    cereal::BinaryInputArchive archive_in(is);
    archive_in(sharded_table);

    cout << "loading reads" << endl;
    loadQueryReads(readFile, readsInfols);

    cout << "processing mapping" << endl;

    string summary;
    summary = matchCidThread(readsInfols, outFile, thread);
    cout << summary << endl;
    return summary;
}

} // namespace CIDmapping