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
using namespace seqan3;

size_t KMERLEN = 6;
size_t BUCKETNUM = 6;
const int K = 25; // K-mer长度
const int BITS_PER_BASE = 2; // 每个碱基使用的比特位数
const int BASES_PER_BYTE = 4; // 每个字节可以存储的碱基数
auto alnconfig = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_score{};

struct Seqinfo {
    bitset<K* BITS_PER_BASE>  encode_seqs;
    uint32_t  xpos;
    uint32_t  ypos;
    template<class Archive>
    void serialize(Archive& ar) {
        ar(encode_seqs, xpos, ypos);
    }
};

struct ReadInfo {
    string id;
    string seq;
    pair<uint64_t, uint64_t> pos;
    string gene;
    string strand;
};

uint32_t murmurHash3(string seq, uint64_t seed) {
    uint32_t out;
    MurmurHash3_x64_128(seq.data(), seq.size(), seed, &out);
    return out;
}

uint32_t min_hash(string seq, size_t len, uint32_t seed, uint8_t& murpos) {
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

std::string reverse_complement(const std::string& sequence) {
    std::string complement(sequence);
    std::reverse(complement.begin(), complement.end());
    for (char& c : complement) {
        switch (c) {
        case 'A':
            c = 'T';
            break;
        case 'T':
            c = 'A';
            break;
        case 'C':
            c = 'G';
            break;
        case 'G':
            c = 'C';
            break;
        default:
            // Handle other characters if necessary
            break;
        }
    }
    return complement;
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

string decode_sequence(const bitset<K* BITS_PER_BASE>& encoded_seq) {
    string decoded_seq;
    for (int i = 0; i < K * BITS_PER_BASE; i += BITS_PER_BASE) {
        int base_bits = (encoded_seq[i] << 1) | encoded_seq[i + 1];
        char base;
        switch (base_bits) {
        case 0: base = 'A'; break;
        case 1: base = 'T'; break;
        case 2: base = 'C'; break;
        case 3: base = 'G'; break;
        }
        decoded_seq += base;
    }
    return decoded_seq;
}

void loadNanoReads(string readFile, vector<ReadInfo>& readsInfols,
    vector<string>  &genels) {
    ifstream inFile(readFile);
    if (!inFile) {
        cerr << "Error: failed to open " <<readFile << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    set<string> gene_unique;
    while (getline(inFile, line)) {
        istringstream ss(line);
        string readid, strand, mpos, mstat, seq, genes, startst, endst;
        getline(ss, readid, '\t');
        getline(ss, strand, '\t');
        getline(ss, mpos, '\t');
        getline(ss, mstat, '\t');
        getline(ss, seq, '\t');
        getline(ss, genes, '\t');
        getline(ss, startst, '\t');
        getline(ss, endst, '\t');
        uint64_t start, end;
        start = stoull(startst);
        end = stoull(endst);
       
        if (genes != "none") {
            stringstream genesS(genes);
            string gene;
            while (getline(genesS, gene, ',')) {
                //cout << gene << endl;
                ReadInfo readInfo;
                readInfo.gene = gene;
                readInfo.id = readid;
                readInfo.seq = seq;
                readInfo.strand = strand;
                readInfo.pos = make_pair(start, end);
                readsInfols.push_back(readInfo);
                if (gene_unique.count(gene) == 0) {
                    //cout << genes << "\t" << gene << endl;
                    gene_unique.insert(gene);
                    genels.push_back(gene);
                }
                
            }
        }
    }
    cout << genels.size() << endl;

}

void loadNanoReads2(string readFile, vector<ReadInfo>& readsInfols,
    vector<string>& genels) {
    ifstream inFile(readFile);
    if (!inFile) {
        cerr << "Error: failed to open " << readFile << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    set<string> gene_unique;
    while (getline(inFile, line)) {
        istringstream ss(line);
        string readid, transcript, seq, genes;
        getline(ss, readid, '\t');
        getline(ss, transcript, '\t');
        getline(ss, seq, '\t');
        getline(ss, genes, '\t');

        if (genes != "none") {
            stringstream genesS(genes);
            string gene;
            while (getline(genesS, gene, ',')) {
                //cout << gene << endl;
                ReadInfo readInfo;
                readInfo.gene = gene;
                readInfo.id = readid;
                readInfo.seq = reverse_complement(seq);
                readInfo.strand = "*";
                readsInfols.push_back(readInfo);
                if (gene_unique.count(gene) == 0) {
                    //cout << genes << "\t" << gene << endl;
                    gene_unique.insert(gene);
                    genels.push_back(gene);
                }

            }
        }
    }
    cout << genels.size() << endl;

}

int jaccardCount(vector<uint64_t> &v1, vector<uint64_t> &v2) {
    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    vector<uint64_t> result;
    set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(result));
    return result.size();
}


int pairwiseAlign(string seq1, string seq2) {
    auto seq1_dna5 = seq1 | seqan3::views::char_to<seqan3::dna5>;
    auto seq2_dna5 = seq2 | seqan3::views::char_to<seqan3::dna5>;
    auto results = seqan3::align_pairwise(std::tie(seq1_dna5, seq2_dna5), alnconfig);
    auto& result = *results.begin();
    return result.score();
}


int editDistance(const std::string& s1, const std::string& s2, int threshold) {
    int m = s1.size(), n = s2.size();
    std::vector<int> dp(n + 1, 0);

    for (int j = 1; j <= n; ++j) dp[j] = j;

    for (int i = 1; i <= m; ++i) {
        int prev = dp[0];
        dp[0] = i;

        for (int j = 1; j <= n; ++j) {
            int temp = dp[j];
            dp[j] = std::min({ dp[j] + 1, dp[j - 1] + 1, prev + (s1[i - 1] != s2[j - 1]) });
            prev = temp;
        }

        if (dp[n] > threshold) return 100; // 如果编辑距离超过阈值，立即返回100
    }

    return dp[n];
}

class ThreadPool {
public:
    ThreadPool(size_t);
    template<class F> void enqueue(F f);
    ~ThreadPool();

private:
    vector<thread> workers;
    queue<function<void()>> tasks;

    mutex queue_mutex;
    condition_variable condition;
    bool stop;
};

inline ThreadPool::ThreadPool(size_t threads)
    : stop(false) {
    for (size_t i = 0; i < threads; ++i)
        workers.emplace_back(
            [this] {
                for (;;) {
                    function<void()> task;

                    {
                        unique_lock<mutex> lock(this->queue_mutex);
                        this->condition.wait(lock,
                            [this] { return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty())
                            return;
                        task = move(this->tasks.front());
                        this->tasks.pop();
                    }

                    task();
                }
            }
    );
}

template<class F>
void ThreadPool::enqueue(F f) {
    {
        unique_lock<mutex> lock(queue_mutex);
        if (stop)
            throw runtime_error("enqueue on stopped ThreadPool");

        tasks.emplace(f);
    }
    condition.notify_one();
}

inline ThreadPool::~ThreadPool() {
    {
        unique_lock<mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for (thread& worker : workers)
        worker.join();
}

void matchCidTask(vector<ReadInfo>& readsInfols, size_t start, size_t end,
    map<string, map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>>& TotalBuckets,
    map<string, map<uint32_t, Seqinfo>>& seqBuckets, ofstream& outputFile, mutex& outputMutex) {

    size_t hitcount = 0;
    for (size_t i = start; i < end; i++) {
        auto startTime = std::chrono::steady_clock::now();
        string seqi = readsInfols[i].seq;
        string strand = readsInfols[i].strand;
        string readid = readsInfols[i].id + "\t" + readsInfols[i].seq;
        string genei = readsInfols[i].gene;
        string resPos, resSeq;

        if (seqi.size() < 20 || seqi.size() > 30) {
            continue;
        }
        uint32_t hashTotal = murmurHash3(seqi, 0);
        auto& bucketi = TotalBuckets[genei];
        auto& bucketseqi = seqBuckets[genei];

        //if exact match;
        //cout << "exact match|" << endl;
        auto it = bucketseqi.find(hashTotal);
        if (it != bucketseqi.end()) {
            Seqinfo& seqinfoi = it->second;
            if (seqinfoi.xpos != 0 && seqinfoi.ypos != 0) {
                //cout << "exactMatch\t" << readid << endl;
                hitcount++;
                resPos = to_string(seqinfoi.xpos) + "_" + to_string(seqinfoi.ypos);
                auto endTime = std::chrono::steady_clock::now();
                std::chrono::duration<double> duration = endTime - startTime;
                lock_guard<mutex> lock(outputMutex);
                outputFile << readid << "\t" << resPos << "\t" << seqi << "\t" << 0 << "\t" << genei << "\t" << strand << "\t" << duration.count() << "\t" << 0 << endl;
            }
        }
        else {
            unordered_set<uint32_t> kmerMap;
            std::unordered_map<uint32_t, int> countMap; // 用于记录满足条件的次数

            for (uint8_t s = 0; s < BUCKETNUM; s++) // 注意这里应该是 s++ 而不是 k++
            {
                uint8_t murpos;
                uint32_t hashMin = min_hash(seqi, KMERLEN, s, murpos);
                std::vector<std::pair<uint32_t, uint8_t>> targetBkt = bucketi[s][hashMin];

                for (size_t tgti = 0; tgti < targetBkt.size(); tgti++)
                {
                    uint32_t seqhashi = targetBkt[tgti].first;
                    uint8_t murposDifference = std::abs(targetBkt[tgti].second - murpos);

                    if (murposDifference < 3)
                    {
                        countMap[seqhashi]++;
                        if (countMap[seqhashi] == 2) // 只有当满足条件两次时才插入到kmerMap
                        {
                            kmerMap.insert(seqhashi);
                        }
                    }
                }
            }
            if (!kmerMap.empty()) {
                int mindistance = 100;
                for (auto it = kmerMap.begin(); it != kmerMap.end(); ++it)
                {
                    uint32_t seqhashi = *it;
                    int distance = 100;
                    string refseqi = decode_sequence(bucketseqi[seqhashi].encode_seqs);
                    //string refseqi = seqHash[hashTotali];
                    EdlibAlignResult matchresult = edlibAlign(seqi.c_str(), seqi.length(), refseqi.c_str(), 25, edlibDefaultAlignConfig());
                    edlibFreeAlignResult(matchresult);
                    if (matchresult.status == EDLIB_STATUS_OK) {
                        distance = matchresult.editDistance;
                    }

                    // cout << "match score: " << mindistance << endl;
                    if (distance < mindistance) {
                        mindistance = distance;
                        resSeq = refseqi;
                        resPos = to_string(bucketseqi[seqhashi].xpos) + "_" + to_string(bucketseqi[seqhashi].ypos);
                    }
                    if (mindistance < 5) {
                        break;
                    }
                }

                if (mindistance != 100) {
                    auto endTime = std::chrono::steady_clock::now();
                    std::chrono::duration<double> duration = endTime - startTime;
                    hitcount++;
                    lock_guard<mutex> lock(outputMutex);
                    outputFile << readid << "\t" << resPos << "\t" << resSeq << "\t" << mindistance << "\t" << genei << "\t" << strand << "\t" << duration.count() << "\t" << kmerMap.size() << endl;

                }
            }
        }
    }
    cout << " hit: " << hitcount << " in " << end - start << endl;
}

void matchCidThread(vector<ReadInfo>& readsInfols,
    map<string, map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>>& TotalBuckets,
    map <string, map<uint32_t, Seqinfo>>& seqBuckets,
    string outputFileName, int numThreads) {

    ThreadPool pool(numThreads);
    mutex outputMutex;
    ofstream outputFile(outputFileName);
    outputFile << "readid\tquerySeq\tcidPos\tquerySeq\teditDi\tgene\tstrand\ttime\tkmsize" << endl;

    size_t n_records = readsInfols.size();
    size_t chunk_size = 100;

    for (size_t start = 0; start < n_records; start += chunk_size) {
        size_t end = min(start + chunk_size, n_records);
        pool.enqueue([&, start, end] {
            matchCidTask(readsInfols, start, end, TotalBuckets, seqBuckets, outputFile, outputMutex);
            });
    }
}

int main(int argc, char* argv[]) {
    string readFile = argv[1];
    string indexFold = argv[2];
    string outFile = argv[3];
    int thread = atoi(argv[4]);
    int errorcount = atoi(argv[5]);

     vector<ReadInfo> readsInfols;
    vector<string> genels;
    vector<string> mygenels;
    map<string, map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>> TotalBuckets;
    map <string, map<uint32_t, Seqinfo>> seqBuckets;

    // load dp index
    std::ifstream is(indexFold, std::ios::binary);
    cereal::BinaryInputArchive archive(is);
    archive(TotalBuckets, seqBuckets);

    cout << "loading reads" << endl;
    loadNanoReads2(readFile, readsInfols, genels);

    cout << "processing mapping" << endl;
    matchCidThread(readsInfols, TotalBuckets, seqBuckets, outFile, thread);

    return 0;
}
