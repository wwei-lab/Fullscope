#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <unordered_map>
#include <mutex>
#include <bitset>
#include <cstdint>
#include <sstream>
#include "edlib.h"
#include <cereal/types/bitset.hpp>  // Include for bitset serialization
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/bitset.hpp>
#include <cereal/types/utility.hpp>  // Add this line for std::pair support
#include <cereal/archives/binary.hpp>

namespace FsUtils {
    struct SeqInfo {
        std::bitset<50> encode_seqs;  // Keep consistent encoding length
        uint32_t xpos;
        uint32_t ypos;

        // Serialization method for cereal
        template<class Archive>
        void serialize(Archive& ar) {
            ar(encode_seqs, xpos, ypos);
        }
    };

// 反向互补转换
std::string reverse_complement(const std::string& seq);

// 字符串分割
std::vector<std::string> split_string(const std::string& s, char delimiter);

// 适配器反向互补预处理
//void precompute_rc(const std::vector<std::pair<std::string, std::string>>& adapters);

// MurmurHash3实现
uint32_t murmur_hash3(const std::string& seq, uint32_t seed);

// Edlib比对包装函数
int edlib_align(const std::string& target, 
                const std::string& query, 
                int max_distance);

std::string join_strings(const std::vector<std::string>& strings, const std::string& delimiter);

// 最小哈希计算
uint32_t min_hash(const std::string& seq, 
                 size_t kmer_len, 
                 uint32_t seed,
                 uint8_t& murpos);

// K-mer哈希生成
std::vector<uint32_t> kmer_hashes(const std::string& seq,
                                  size_t kmer_len,
                                  uint32_t seed);

// 序列编码（25-mer使用50位存储）
std::bitset<50> encode_sequence(const std::string& seq);

// 序列解码
std::string decode_sequence(const std::bitset<50>& encoded);

} // namespace FsUtils

#endif // UTILS_H