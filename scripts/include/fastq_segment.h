#ifndef FASTQ_PROCESSOR_H
#define FASTQ_PROCESSOR_H

#include <vector>
#include <string>
#include <utility>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <atomic>
#include <thread>
#include "utils.h"
using namespace std;
namespace FsFastqSegment {

// 比对结果数据结构
struct AlignmentResult {
    int start;
    int end;
    int distance;
    std::string adapter_id;
    bool is_rc;
};

// 分段结果数据结构
struct Segment {
    int start;
    int end;
    std::string adapter_id;
    std::string strand;
    std::string read_id;
};

// 读取数据结构
struct ReadData {
    std::string header;
    std::string sequence;
    std::string quality;
    size_t read_number;
};

/**
 * @brief 在目标序列中查找适配器
 * @param target 目标序列
 * @param adapters 适配器列表（ID，序列）
 * @param threshold 匹配阈值（错误率）
 * @return 比对结果列表
 */
std::vector<AlignmentResult> find_adapters(
    const std::string& target,
    const std::vector<std::pair<std::string, std::string>>& adapters,
    double threshold);

/**
 * @brief 根据比对结果生成分段信息
 * @param matches 比对结果列表
 * @param header 读段头信息
 * @param seq_length 序列总长度
 * @return 分段结果列表
 */
std::vector<Segment> generate_segments(
    const std::vector<AlignmentResult>& matches,
    const std::string& header,
    int seq_length);

/**
 * @brief 多线程处理FASTQ文件
 * @param fastq_path 输入文件路径
 * @param adapters 适配器列表（ID，序列）
 * @param threshold 匹配阈值（错误率）
 * @param thread_count 线程数量
 * @param output_path 输出文件路径
 */

std::vector<std::pair<std::string, std::string>> load_adapters(const std::string& path);
void precompute_rc(const std::vector<std::pair<std::string, std::string>>& adapters);
string fastq_segment_main(const string& fastq_path,
    const string& adapter_path,
    double threshold,
    unsigned int thread_count,
    const string& output_path);
} // namespace FastqProcessor

#endif // FASTQ_PROCESSOR_H