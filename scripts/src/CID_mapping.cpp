#include "CID_mapping.h"
#include <chrono>
#include <thread>
#include <unordered_set>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>

using namespace std;
namespace CIDmapping {

size_t KMERLEN = 6;
size_t BUCKETNUM = 6;
const int K = 25;
const int BITS_PER_BASE = 2;
const int BASES_PER_BYTE = 4;
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
        string readid, strand, pos, status, seq, genes;
        getline(ss, readid, '\t');
        getline(ss, strand, '\t');
        getline(ss, pos, '\t');
        getline(ss, status, '\t');
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
                readInfo.seq = seq;
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

int pairwiseAlign(string seq1, string seq2) {
    auto seq1_dna5 = seq1 | seqan3::views::char_to<seqan3::dna5>;
    auto seq2_dna5 = seq2 | seqan3::views::char_to<seqan3::dna5>;
    auto results = seqan3::align_pairwise(std::tie(seq1_dna5, seq2_dna5), alnconfig);
    auto& result = *results.begin();
    return result.score();
}

void matchCidChunk(vector<ReadInfo>& readsInfols,
    map<string, map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>>& TotalBuckets,
    map<string, map<uint32_t, Seqinfo>>& seqBuckets, 
    unordered_map<string, int>& gene_threshold_map,
    ofstream& outputFile, int start, int end, 
    mutex& mtx) {
    size_t hitcount = 0;
    vector<string> outputBuffer;
    const size_t bufferSize = 10000; 

    //cout << "processing" << endl;
    for (int i = start; i <= end; i++) {
        auto startTime = std::chrono::steady_clock::now();
        string seqi = readsInfols[i].seq;
        string strand = readsInfols[i].strand;
        string readid = readsInfols[i].id + "\t" + readsInfols[i].seq;
        string genei = readsInfols[i].gene;
        string resPos, resSeq;
        auto it_ref = gene_threshold_map.find(genei);
        int ref_threshold = 0;
        if (it_ref != gene_threshold_map.end()) {
            ref_threshold = it_ref->second;
        }

        if (seqi.size() < 20 || seqi.size() > 30) {
            continue;
        }
        uint32_t hashTotal = murmurHash3(seqi, 0);
        auto& bucketi = TotalBuckets[genei];
        auto& bucketseqi = seqBuckets[genei];

        //cout << "hash match" << endl;
        // Exact match
        auto it = bucketseqi.find(hashTotal);
        if (it != bucketseqi.end()) {
            Seqinfo& seqinfoi = it->second;
            if (seqinfoi.xpos != 0 && seqinfoi.ypos != 0) {
                hitcount++;
                resPos = to_string(seqinfoi.xpos) + "_" + to_string(seqinfoi.ypos);
                auto endTime = std::chrono::steady_clock::now();
                std::chrono::duration<double> duration = endTime - startTime;
                outputBuffer.emplace_back(readid + "\t" + resPos + "\t" + seqi + "\t0\t" + to_string(ref_threshold) + "\t" + 
                genei + "\t" + strand + "\t" + to_string(duration.count()) + "\t0\t0");

                continue;
            }
        }

        // kmer match
        unordered_set<uint32_t> kmerMap;
        std::unordered_map<uint32_t, int> countMap;
        //cout << "kmer match" << endl;
        for (uint8_t s = 0; s < BUCKETNUM; s++) {
            uint8_t murpos;
            uint32_t hashMin = min_hash(seqi, KMERLEN, s, murpos);
            std::vector<std::pair<uint32_t, uint8_t>> targetBkt = bucketi[s][hashMin];

            for (size_t tgti = 0; tgti < targetBkt.size(); tgti++) {
                uint32_t seqhashi = targetBkt[tgti].first;
                uint8_t murposDifference = std::abs(targetBkt[tgti].second - murpos);

                if (murposDifference < 3) {
                    countMap[seqhashi]++;
                    if (countMap[seqhashi] == 2) {
                        kmerMap.insert(seqhashi);
                    }
                }
            }
        }

        auto endTimekmer = std::chrono::steady_clock::now();
        std::chrono::duration<double> durationkmer = endTimekmer - startTime;

        if (!kmerMap.empty()) {
            int mindistance = 100;
            for (auto it = kmerMap.begin(); it != kmerMap.end(); ++it) {
                uint32_t seqhashi = *it;
                int distance = 100;
                string refseqi = decode_sequence(bucketseqi[seqhashi].encode_seqs);
                EdlibAlignResult matchresult = edlibAlign(seqi.c_str(), seqi.length(), refseqi.c_str(), 25, edlibDefaultAlignConfig());
                edlibFreeAlignResult(matchresult);
                if (matchresult.status == EDLIB_STATUS_OK) {
                    distance = matchresult.editDistance;
                }

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
                outputBuffer.emplace_back(readid + "\t" + resPos + "\t" + resSeq + "\t" + to_string(mindistance) + "\t" + to_string(ref_threshold) + 
                "\t" + genei + "\t" + strand + "\t" + to_string(duration.count()) + "\t" + to_string(durationkmer.count()) + "\t" + to_string(kmerMap.size()));
            }
            
        }

        if (outputBuffer.size() >= bufferSize) {
            mtx.lock();
            for (const auto& line : outputBuffer) {
                outputFile << line << endl;
            }
            outputBuffer.clear();
            mtx.unlock();
        }
    }

    mtx.lock();
    for (const auto& line : outputBuffer) {
        outputFile << line << endl;
    }
    mtx.unlock();

    cout << "hit: " << hitcount << " in " << end - start << endl;
}

void matchCidThread(vector<ReadInfo>& readsInfols,
    map<string, map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>>& TotalBuckets,
    map <string, map<uint32_t, Seqinfo>>& seqBuckets,
    unordered_map<string, int>& gene_threshold_map,
    string outputFileName, int numThreads) {
    vector<thread> threads;
    mutex mtx;
    int n_records = readsInfols.size();
    int chunk_size = n_records / numThreads;
    int remainder = n_records % numThreads;

    ofstream outputFile(outputFileName);
    outputFile << "readid\tquerySeq\tcidPos\trefSeq\teditDi\teditDi_ref\tgene\tstrand\ttime\tkmTime\tkmSize" << endl;
    int start = 0, end = 0;
    for (int i = 0; i < numThreads; i++) {
        start = end;
        end = start + chunk_size;
        if (i == numThreads - 1) {
            end += remainder - 1;
        }
        cout << start << "-" << end << endl;
        threads.emplace_back(matchCidChunk, ref(readsInfols), ref(TotalBuckets), ref(seqBuckets),  ref(gene_threshold_map),
        ref(outputFile), start, end, ref(mtx));
    }

    for (auto& t : threads) {
        t.join();
    }
}

string filterResults(const string& inputFile, const string& outputFile) {
    ifstream in(inputFile);
    ofstream out(outputFile);
    string line;
    
    // 读取并写入标题行
    getline(in, line);
    out << line << "\n";

    int total = 0;
    int passed = 0;
    vector<string> buffer;
    const size_t bufferSize = 10000;

    while (getline(in, line)) {
        vector<string> fields;
        string field;
        istringstream iss(line);
        
        // 解析TSV字段
        while (getline(iss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 10) {
            cerr << "Invalid line format: " << line << endl;
            continue;
        }

        int editDi, editDi_ref;
        try {
            editDi = stoi(fields[4]);
            editDi_ref = stoi(fields[5]);
        } catch (...) {
            cerr << "Error parsing numbers in line: " << line << endl;
            continue;
        }

        total++;

        // 过滤条件
        if (editDi <= editDi_ref) {
            passed++;
            buffer.push_back(line);
            
            // 缓冲写入
            if (buffer.size() >= bufferSize) {
                for (const auto& l : buffer) out << l << "\n";
                buffer.clear();
            }
        }
    }

    // 写入剩余内容
    if (!buffer.empty()) {
        for (const auto& l : buffer) out << l << "\n";
    }

    // 计算并输出统计信息
    double ratio = 0;
    if (total > 0) {
        ratio = static_cast<double>(passed) * 100 / total;
        cout << "Passing percentage: " << ratio << "%" << endl;
    } else {
        cout << "No sequences processed" << endl;
    }
    string summary;
    summary = "Total sequences: " + to_string(total) + "\t" + "Passed sequences: " + to_string(passed) + "Passed ratio: " + to_string(ratio) + "\n";

    in.close();
    out.close();
    return summary;
}

string CIDmap_processer(string readFile, string indexFold, int thread, int kmerlen, int bucketnum, string outFile, string thresholdFile) {
    KMERLEN = kmerlen;
    BUCKETNUM = bucketnum;

    vector<ReadInfo> readsInfols;
    vector<string> genels;
    vector<string> mygenels;
    map <string, map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>> TotalBuckets;
    map <string, map<uint32_t, Seqinfo>> seqBuckets;

    // load dp index
    cout << "loading index" << endl;
    std::ifstream is(indexFold, std::ios::binary);
    cereal::BinaryInputArchive archive(is);
    archive(TotalBuckets, seqBuckets);
     // load dynamic threshold
    unordered_map<string, int> gene_threshold_map;
    ifstream thresholdStream(thresholdFile);
    string line;
    getline(thresholdStream, line); // skip column
    while (getline(thresholdStream, line)) {
        istringstream iss(line);
        string gene;
        int count, threshold;
        if (iss >> gene >> count >> threshold) {
            gene_threshold_map[gene] = threshold;
        }
    }

    cout << "loading reads" << endl;
    loadNanoReads(readFile, readsInfols, genels);

    cout << "processing mapping" << endl;
    matchCidThread(readsInfols, TotalBuckets, seqBuckets, gene_threshold_map, outFile, thread);

    string summary;
    cout << "\nStart filtering results..." << endl;
    summary = filterResults(outFile, outFile + ".filtered");
    cout << "Filtering completed.\n";
    return summary;
}

} // namespace CIDmapping