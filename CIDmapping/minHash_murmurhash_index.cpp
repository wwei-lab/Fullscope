#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include "MurmurHash3.h"
#include<time.h>
#include<map>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/bitset.hpp>  // for std::bitset
#include <htslib/sam.h>
#include <set>
#include <thread>
#include <mutex>
#include <sstream>
#include <unistd.h> 
#include <bitset>
using namespace std;

size_t KMERLEN = 6;
size_t BUCKETNUM = 6;
const int K = 25; // K-mer长度
const int BITS_PER_BASE = 2; // 每个碱基使用的比特位数
const int BASES_PER_BYTE = 4; // 每个字节可以存储的碱基数

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
    iss >> pages; // 获取内存页数

    // 内存页大小（以字节为单位）
    const std::size_t page_size = sysconf(_SC_PAGESIZE);
    return pages * page_size; // 内存使用量
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

uint32_t murmurHash3(string seq, uint32_t seed) {
    uint32_t out;
    MurmurHash3_x64_128(seq.data(), seq.size(), seed, &out);
    return out;
}

//get min hash of different k-mer
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

void kmer_hash(string seq, size_t len, uint32_t seed, vector<uint32_t>& hashvec) {
    for (size_t i = 0; i <= (seq.size() - len); i++)
    {
        string seqsub = seq.substr(i, len);
        uint32_t hasht = murmurHash3(seqsub, seed);
        hashvec.push_back(hasht);
    }
}

// build buckets for sequence vector
void bucket_build(map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>& bucketi, string seq, uint32_t hashTotal) {
    uint8_t murpos;
    for (uint8_t s = 0; s < BUCKETNUM; s++)
    {
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


void map_stat(map<uint32_t, vector<uint32_t>> &myHash) {
    int numKeys = myHash.size();
    std::vector<int> value_sizes;
    for (auto const& [key, value] : myHash) {
        value_sizes.push_back(value.size());
    }
    std::sort(value_sizes.begin(), value_sizes.end());
    int min_size = value_sizes.front();
    int median_size = value_sizes[value_sizes.size() / 2];
    int max_size = value_sizes.back();
    std::cout << "keys: " << numKeys <<  "\tMin: " << min_size << "\tMedian: " << median_size << "\tMax: " << max_size << std::endl;
}

void gene_index_precess_chunk(vector<pair <string, string>>& geneRanges, string bamFile, 
    uint64_t gstart, uint64_t gend, 
    map<string, map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>>& TotalBuckets,
    map <string, map<uint32_t, Seqinfo>> &seqBuckets,
    mutex & mtx) {
    string cidseq, cidpos;
    string rqns;
    uint32_t hashTotal;
    // extract cid seq and pos
    set<uint32_t> uniqueSet;
    // extract read align pos
    //vector<map<size_t, map<uint32_t, vector<uint32_t>>>> buckets;
    map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>> bucketi;
    map<uint32_t, Seqinfo> seqseti;
    Seqinfo seqinfoi;
    samFile* fp = sam_open(bamFile.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(fp);
    hts_idx_t* idx = sam_index_load(fp, bamFile.c_str());
    bam1_t* record = bam_init1();

    for (int gi = gstart; gi <= gend; gi++) {
        for (auto& outerPair : bucketi) {
            // Iterate through the inner map
            for (auto& innerPair : outerPair.second) {
                // Clear and shrink the vector
                innerPair.second.clear();
                innerPair.second.shrink_to_fit();
            }
            // Clear the inner map
            outerPair.second.clear();
        }
        // Clear the outer map
        bucketi.clear();
        seqseti.clear();

        //buckets.clear();
        uniqueSet.clear();
        //cout << "extract gene" << endl;
        string geneid = geneRanges[gi].first;
        string rg = geneRanges[gi].second;
        const char* c_rg = rg.c_str();
        hts_itr_t* iter = sam_itr_querys(idx, header, c_rg);
        
        //cout << "process reads" << endl;
        while (sam_itr_next(fp, iter, record) >= 0) {
            rqns = bam_get_qname(record);
            cidseq = rqns.substr(37, 25);
            for (int k = 0; rqns[37 + 28 + k] != '|'; k++) {
                cidpos = cidpos + rqns[37 + 28 + k];
            }
            //cout << rqns << endl;
            hashTotal = murmurHash3(cidseq, 0);
            //cout << hashTotal << endl;
            if (uniqueSet.count(hashTotal) == 0) {
                uniqueSet.insert(hashTotal);
                bucket_build(bucketi, cidseq, hashTotal);
                seqinfoi.encode_seqs = encode_sequence(cidseq);

                size_t pos = cidpos.find('_');
                uint32_t xpos = std::stoi(cidpos.substr(5, pos));
                uint32_t ypos = std::stoi(cidpos.substr(pos + 1));
                seqinfoi.xpos = xpos;
                seqinfoi.ypos = ypos;
                seqseti[hashTotal] = seqinfoi;
                //mtx.lock();
                //write to read and cid content file
                //outputFile << geneid << "\t" << cidseq << "\t" << cidpos << "\t" << hashTotal << endl;
                //mtx.unlock();
                //cout << "build bucket" << endl;
                /*
                *                  for (size_t seed = 0; seed < bucketNum; seed++)
                {
                    
                    buckets.push_back(bucketi);
                }
                */
            }
            cidseq = "";
            cidpos = "";
        }
        hts_itr_destroy(iter);

        if (uniqueSet.size() == 0) {
            continue;
        }
        TotalBuckets[geneid] = bucketi;
        seqBuckets[geneid] = seqseti;
        std::cout << geneid << " unique_CIDseq_size\t" << uniqueSet.size() << std::endl;

        //kmerHash.clear();
        //rangeHash.clear();
        std::size_t memory_usage = getCurrentMemoryUsage();
        std::cout << "Current memory usage: " << memory_usage / (1024 * 1024) << " MB" << std::endl;

    }
    bam_destroy1(record);
    hts_idx_destroy(idx);
    sam_close(fp);
    bam_hdr_destroy(header);

}

void gene_index_process(vector<pair <string, string>>& geneRanges, string bamFile, int numThreads,
    map<string, map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>>& TotalBuckets,
    map <string, map<uint32_t, Seqinfo>>& seqBuckets) {
    mutex mtx;
    uint32_t genecount = geneRanges.size();
    uint32_t chunk_size = genecount / numThreads;
    uint32_t remainder = genecount % numThreads;
    vector<thread> threads;
    uint64_t start = 0, end = 0;
    for (uint64_t i = 0; i < numThreads; i++) {
        start = end;
        end = start + chunk_size;
        if (i == numThreads - 1) {
            end += remainder;
        }
        cout << start << "-" << end << endl;
        threads.emplace_back(gene_index_precess_chunk, ref(geneRanges), bamFile, start, end - 1,ref(TotalBuckets), ref(seqBuckets), ref(mtx));
    }

    for (auto& t : threads) {
        t.join();
    }

}

int main(int argc, char* argv[]) {
    string bamFile = argv[1];
    string geneFile = argv[2];
    string outFold = argv[3];
    int numThreads = atoi(argv[4]);
    KMERLEN = atoi(argv[5]);
/*
*    if (!filesystem::exists(outFold)) // if dir exist
    {
        if (filesystem::create_directory(outFold)) // create dir
        {
            cout << "Create fold" << std::endl;
        }
        else
        {
            cerr << "Create fold fail！" << std::endl;
        }
    }
*/
 

    map<string, map<uint8_t, map<uint32_t, vector<pair<uint32_t, uint8_t>>>>> TotalBuckets;
    map <string, map<uint32_t, Seqinfo>> seqBuckets;
    //map<uint64_t, string> seqHash;
    vector<pair <string, string>> geneRanges;
    cout << "loading readlist" << endl;
    loadGeneRanges(geneFile, geneRanges);
    cout << "building index" << endl;

    gene_index_process(geneRanges, bamFile, numThreads, TotalBuckets, seqBuckets);

    //output index
    string outFm = outFold +  ".bin";
    std::ofstream os{ outFm, std::ios::binary };
    cereal::BinaryOutputArchive oarchive{ os };
    oarchive(TotalBuckets, seqBuckets);
    os.close();

    cout << "All index has been built successfully !" << endl;
    return 0;
}
