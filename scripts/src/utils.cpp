#include "utils.h"
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <cctype>
#include "edlib.h"
#include "MurmurHash3.h"
using namespace std;

namespace FsUtils {
    // 反向互补转换
    string reverse_complement(const string& seq) {
        static const unordered_map<char, char> rc_map {
            {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'},
            {'a','t'}, {'t','a'}, {'c','g'}, {'g','c'},
            {'N','N'}, {'n','n'}
        };
        string rc(seq.rbegin(), seq.rend());
        for(auto& c : rc) {
            c = rc_map.count(c) ? rc_map.at(c) : 'N';
        }
        return rc;
    }

    // 字符串分割
    vector<string> split_string(const string& s, char delimiter) {
        vector<string> tokens;
        string token;
        //istringstream tokenStream(s);
        std::istringstream tokenStream(s);
        while (getline(tokenStream, token, delimiter)) {
            if (!token.empty()) {
                tokens.push_back(token);
            }
        }
        return tokens;
    }

    std::string join_strings(const std::vector<std::string>& strings, const std::string& delimiter) {
        std::ostringstream oss;
        for (size_t i = 0; i < strings.size(); ++i) {
            if (i != 0) oss << delimiter;
            oss << strings[i];
        }
        return oss.str();
    }
    // 适配器反向互补预处理
    // void precompute_rc(const vector<pair<string, string>>& adapters) {
    //     static unordered_map<string, string> rc_cache;
    //     static mutex cache_mutex;
        
    //     lock_guard<mutex> lock(cache_mutex);
    //     for(const auto& adapter : adapters) {
    //         const string& seq = adapter.second;
    //         if(rc_cache.count(seq) == 0) {
    //             rc_cache[seq] = reverse_complement(seq);
    //         }
    //     }
    // }

    // MurmurHash3实现
    uint32_t murmur_hash3(const string& seq, uint32_t seed) {
        uint32_t out;
        MurmurHash3_x64_128(seq.data(), seq.size(), seed, &out);
        return out;
    }

    // Edlib比对包装函数
    int edlib_align(const string& target, const string& query, int max_distance) {
        //EdlibAlignConfig config = edlibNewAlignConfig(max_distance, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE);
        EdlibAlignResult result = edlibAlign(query.c_str(), query.length(),
                                            target.c_str(), target.length(),
                                            edlibDefaultAlignConfig());
        int distance = (result.status == EDLIB_STATUS_OK) ? result.editDistance : -1;
        edlibFreeAlignResult(result);
        return distance;
    }

    // 最小哈希计算
    uint32_t min_hash(const string& seq, size_t kmer_len, uint32_t seed, uint8_t& murpos) {
        if(seq.length() < kmer_len) return UINT32_MAX;
        
        uint32_t min_val = UINT32_MAX;
        murpos = 0;

        for(uint8_t i = 0; i <= seq.size() - kmer_len; ++i) {
            string kmer = seq.substr(i, kmer_len);
            uint32_t hash_val = murmur_hash3(kmer, seed);
            
            if(hash_val < min_val) {
                min_val = hash_val;
                murpos = i;
            }
        }
        return min_val;
    }

    // K-mer哈希生成
    vector<uint32_t> kmer_hashes(const string& seq, size_t kmer_len, uint32_t seed) {
        vector<uint32_t> hashes;
        if(seq.length() < kmer_len) return hashes;

        for(size_t i = 0; i <= seq.size() - kmer_len; ++i) {
            string kmer = seq.substr(i, kmer_len);
            hashes.push_back(murmur_hash3(kmer, seed));
        }
        return hashes;
    }

}