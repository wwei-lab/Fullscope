#include "fastq_segment.h"
#include "utils.h"
#include <algorithm>
#include <thread>
#include <mutex>
#include <queue>
#include <fstream>
#include <sstream>
#include <numeric>
#include <vector>
#include <atomic>  
using namespace std;
using namespace FsUtils;

namespace FsFastqSegment {

// 线程控制参数
static unsigned int THREAD_COUNT = thread::hardware_concurrency();
unordered_map<string, string> rc_cache;

vector<AlignmentResult> find_adapters(const string& target,
                                    const vector<pair<string, string>>& adapters,
                                    double threshold) 
{
    // 建立k-mer索引 (k=10)
    const int k = 10;
    unordered_map<string, vector<int>> kmer_index;
    if (target.length() >= k) {
        for (int i = 0; i <= static_cast<int>(target.length()) - k; ++i) {
            string kmer = target.substr(i, k);
            kmer_index[kmer].push_back(i);
        }
    }

    vector<AlignmentResult> results;
    unordered_set<int> processed_positions;  // 用于记录已处理的位置

    for (const auto& adapter : adapters) {
        const string& id = adapter.first;
        const string& seq = adapter.second;
        const int L = seq.length();
        const int max_ed = static_cast<int>(threshold * L);
        
        // 如果序列太短或阈值太高，退化为暴力搜索
        if (L < k || max_ed >= L - 1) {
            for (int pos = 0; pos <= static_cast<int>(target.length()) - L; ++pos) {
                string window = target.substr(pos, L);
                
                // 正向比对
                int dist = edlib_align(window, seq, max_ed);
                if (dist != -1 && dist <= max_ed) {
                    results.push_back({pos, pos + L - 1, dist, id, false});
                }
                
                // 反向比对
                int rc_dist = edlib_align(window, rc_cache[seq], max_ed);
                if (rc_dist != -1 && rc_dist <= max_ed) {
                    results.push_back({pos, pos + L - 1, rc_dist, id + "_RC", true});
                }
            }
        } 
        else {
            // 使用k-mer索引优化
            unordered_set<int> candidate_positions;
            
            // 查找候选位置函数
            auto find_candidates = [&](const string& sequence) {
                const int step = max_ed + 1;  // 鸽巢原理步长
                for (int i = 0; i <= L - k; i += step) {
                    string kmer = sequence.substr(i, k);
                    auto it = kmer_index.find(kmer);
                    if (it != kmer_index.end()) {
                        for (int pos : it->second) {
                            int start_pos = pos - i;
                            if (start_pos >= 0 && start_pos <= static_cast<int>(target.length()) - L) {
                                candidate_positions.insert(start_pos);
                            }
                        }
                    }
                }
            };
            
            // 查找正向序列候选位置
            find_candidates(seq);
            
            // 查找反向互补序列候选位置
            find_candidates(rc_cache[seq]);
            
            // 处理候选位置
            for (int pos : candidate_positions) {
                // 跳过已处理位置（不同adapter可能在同一位置）
                if (processed_positions.find(pos) != processed_positions.end()) 
                   continue;
                
                string window = target.substr(pos, L);
                
                // 正向比对
                int dist = edlib_align(window, seq, max_ed);
                if (dist != -1 && dist <= max_ed) {
                    results.push_back({pos, pos + L - 1, dist, id, false});
                    processed_positions.insert(pos);
                }
                
                // 反向比对
                int rc_dist = edlib_align(window, rc_cache[seq], max_ed);
                if (rc_dist != -1 && rc_dist <= max_ed) {
                    results.push_back({pos, pos + L - 1, rc_dist, id + "_RC", true});
                    processed_positions.insert(pos);
                }
            }
        }
    }

    // 结果过滤 (保持原逻辑)
    sort(results.begin(), results.end(), [](const AlignmentResult& a, const AlignmentResult& b) {
        return a.start != b.start ? a.start < b.start : a.distance < b.distance;
    });

    vector<AlignmentResult> filtered;
    for (const auto& res : results) {
        if (!filtered.empty() && res.start <= filtered.back().end) {
            if (res.distance < filtered.back().distance) {
                filtered.pop_back();
                filtered.push_back(res);
            }
        } else {
            filtered.push_back(res);
        }
    }
    return filtered;
}

std::vector<Segment> generate_segments(
    const std::vector<AlignmentResult>& matches,
    const std::string& header,
    int seq_length)
{
    std::vector<Segment> segments;
    for(size_t j = 0; j < matches.size(); ++j) {
        const auto& curr = matches[j];
        if(!curr.is_rc) {
            if(j == 0) {
                int seg_start = 1;
                int seg_end = curr.start + 1;
                if(seg_start < seg_end) {
                    segments.push_back({seg_start, seg_end, curr.adapter_id, "+", header});
                }
            }
            if(j > 0) {
                const auto& prev = matches[j-1];
                if(!prev.is_rc) {
                    int seg_start = prev.end + 1;
                    int seg_end = curr.start + 1;
                    if(seg_start < seg_end) {
                        string combined = prev.adapter_id + "-" + curr.adapter_id;
                        segments.push_back({seg_start, seg_end, combined, "+", header});
                    }
                }
            }
        } else {
            //cout << "strand -" << endl;
            if(j == matches.size() - 1) {
                int seg_start = curr.end + 1;
                int seg_end = seq_length;
                if(seg_start <= seg_end) {
                    segments.push_back({seg_start, seg_end, curr.adapter_id, "-", header});
                }
            } else {
                
                const auto& next = matches[j+1];
                if(next.is_rc) {
                    int seg_start = curr.end + 1;
                    int seg_end = next.start + 1;
                    if(seg_start < seg_end) {
                        string combined = next.adapter_id + "-" + curr.adapter_id;
                        segments.push_back({seg_start, seg_end, combined, "-", header});
                    }
                }
            }
        }
    }
    return segments;
}

void precompute_rc(const vector<pair<string, string>>& adapters) {
    //unordered_map<string, string> rc_cache;
    static unordered_map<char, char> rc_map {
        {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'},
        {'a','t'}, {'t','a'}, {'c','g'}, {'g','c'},
        {'N','N'}, {'n','n'}
    };
    for(const auto& adapter : adapters) {
        const string& seq = adapter.second;
        string rc(seq.rbegin(), seq.rend());
        for(auto& c : rc) c = rc_map.count(c) ? rc_map[c] : 'N';
        rc_cache[seq] = rc;
    }
}

vector<pair<string, string>> load_adapters(const string& path) {
    vector<pair<string, string>> adapters;
    ifstream file(path);
    if(!file.is_open()) {
        cerr << "Error opening adapter file: " << path << endl;
        exit(1);
    }

    string line, current_id, seq;
    while(getline(file, line)) {
        if(line.empty()) continue;
        if(line[0] == '>') {
            if(!current_id.empty()) {
                adapters.emplace_back(current_id, seq);
                seq.clear();
            }
            size_t space_pos = line.find(' ');
            line.erase(remove(line.begin(), line.end(), '\r'), line.end());
            current_id = (space_pos != string::npos) ? line.substr(1, space_pos-1) : line.substr(1);
        } else {
            line.erase(remove(line.begin(), line.end(), '\r'), line.end());
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            seq += line;
        }
    }
    if(!current_id.empty()) {
        adapters.emplace_back(current_id, seq);
    }
    return adapters;
}


// 处理单个任务块的函数
void handle_fastq_chunk(
    string& fastq_path,
    vector<pair<string, string>>& adapters,
    double threshold,
    uint64_t start_idx,
    uint64_t end_idx,
    string& tmp_output_path,
    string& tmp_summary_path,
    size_t& local_segments,
    vector<streampos>& read_offsets
) {
    //cout <<  " processing "
    //    << start_idx << "-" <<  end_idx << " reads\n";
    ofstream out_file(tmp_output_path, ios::binary);
    if (!out_file) {
        cerr << "Error opening: " << tmp_output_path << endl;
        return;
    }
    ofstream summary_file(tmp_summary_path, ios::binary);
    // char out_buffer[1024 * 1024];  // 1MB输出缓冲区
    // out_file.rdbuf()->pubsetbuf(out_buffer, sizeof(out_buffer));
    
    // 打开输入文件并定位
    ifstream in_file(fastq_path);
    // char in_buffer[1024 * 1024];  // 1MB输入缓冲区
    // in_file.rdbuf()->pubsetbuf(in_buffer, sizeof(in_buffer));
    
    // 定位到任务块起始位置
    in_file.seekg(read_offsets[start_idx]);
    
    // 处理当前任务块
    for (uint64_t current_read = start_idx; current_read < end_idx; ++current_read) {
        string header, sequence, plus, quality;
        if (!getline(in_file, header)) break;
        getline(in_file, sequence);
        getline(in_file, plus);
        getline(in_file, quality);
        
        // 生成全局唯一ID
        string global_id = "read_" + to_string(current_read);
        
        // 比对和分段
        auto matches = find_adapters(sequence, adapters, threshold);
        if (matches.empty()) {
            out_file << "@" << global_id << "|Noadapter\n"
                     << sequence << "\n+\n" << quality << "\n";
            
            summary_file << global_id << "\tNoadapter\t0\t"
                         << sequence.size() << "\t*\n";
            local_segments++;
        } else {
            auto segments = generate_segments(matches, global_id, sequence.length());
            for (const auto& seg : segments) {
                int start_idx = seg.start - 1;
                int length = seg.end - seg.start + 1;
                if (length < 100) continue;
                if (start_idx < 0 || start_idx + length > sequence.size()) continue;
                
                string subseq, subqual;
                if (seg.strand == "+") {
                    int reslen = sequence.size() - start_idx - length;
                    int extendlen = min(20, reslen);
                    subseq = sequence.substr(start_idx, length + extendlen);
                    subqual = quality.substr(start_idx, length + extendlen);
                } else {
                    int extendlen = min(20, start_idx);
                    subseq = sequence.substr(start_idx - extendlen, length + extendlen);
                    subqual = quality.substr(start_idx - extendlen, length + extendlen);
                    subseq = reverse_complement(subseq);
                    reverse(subqual.begin(), subqual.end());
                }
                
                out_file << "@" << seg.read_id << "|" << seg.adapter_id << "|"
                         << seg.start << "-" << seg.end << "(" << seg.strand << ")\n"
                         << subseq << "\n+\n" << subqual << "\n";
                
                summary_file << seg.read_id << "\t" << seg.adapter_id << "\t"
                             << seg.start << "\t" << seg.end << "\t" << seg.strand << "\n";
                local_segments++;
            }
        }
    }
    summary_file.close();
    out_file.close();
}

// 主控函数
string process_fastq(string fastq_path,
                     vector<pair<string, string>>& adapters,
                     double threshold,
                     unsigned int thread_count,
                     string output_path) 
{
    // ===== 1. 预计算文件偏移量 =====
    vector<streampos> read_offsets;
    {
        ifstream in_file(fastq_path);
        string line;
        streampos pos = in_file.tellg();
        
        while (getline(in_file, line)) {
            if (!line.empty() && line[0] == '@') {
                read_offsets.push_back(pos);
                // 跳过序列和质量行
                for (int i = 0; i < 3; ++i) getline(in_file, line);
            }
            pos = in_file.tellg();
        }
    }
    size_t total_reads = read_offsets.size();

    // ===== 2. 任务分块 =====
    vector<pair<size_t, size_t>> chunks;
    const size_t chunk_size = total_reads / thread_count;
    size_t remaining = total_reads % thread_count;
    size_t start_idx = 0;
    
    for (size_t i = 0; i < thread_count; ++i) {
        size_t end_idx = start_idx + chunk_size + (i < remaining ? 1 : 0);
        chunks.emplace_back(start_idx, end_idx);
        start_idx = end_idx;
    }

    // ===== 3. 多线程处理 =====
    vector<thread> workers;
    vector<string> tmp_outputs;
    vector<string> tmp_summaries;
    vector<size_t> segment_counts(thread_count, 0);
    size_t total_segments = 0;
    //atomic<size_t> total_segments(0);

    // 预生成所有临时文件路径
    for (size_t tid = 0; tid < thread_count; ++tid) {
        tmp_outputs.push_back(output_path + ".part" + to_string(tid) + ".fq");
        tmp_summaries.push_back(output_path + ".part" + to_string(tid) + ".summary.tsv");
    }

    for (size_t tid = 0; tid < thread_count; ++tid) {
        workers.emplace_back([&, tid] {
            size_t local_segments = 0;  // 线程局部变量
            
            handle_fastq_chunk(
                fastq_path,
                adapters,
                threshold,
                chunks[tid].first,
                chunks[tid].second,
                tmp_outputs[tid],
                tmp_summaries[tid],
                local_segments,
                read_offsets
            );
            
            segment_counts[tid] = local_segments;
        });
    }

    // ===== 4. 等待线程完成 =====
    for (auto& t : workers) {
        t.join();
    }

    // ===== 5. 高效合并文件 =====
    ofstream final_out(output_path);
    ofstream final_summary(output_path + ".summary.tsv");
    // char merge_buffer[1024 * 1024];  // 1MB合并缓冲区
    // final_out.rdbuf()->pubsetbuf(merge_buffer, sizeof(merge_buffer));
    
    for (size_t i = 0; i < thread_count; ++i) {
        total_segments = total_segments + segment_counts[i];
        // 合并主输出
        ifstream part_out(tmp_outputs[i]);
        final_out << part_out.rdbuf();
        part_out.close();
        remove(tmp_outputs[i].c_str());
        
        // 合并摘要
        ifstream part_summary(tmp_summaries[i]);
        final_summary << part_summary.rdbuf();
        part_summary.close();
        remove(tmp_summaries[i].c_str());
    }
    final_out.close();
    final_summary.close();
    // ===== 6. 生成统计摘要 =====
    return "Total reads: " + to_string(total_reads) + 
           "\tProcessed segments: " + to_string(total_segments) + "\n";
}

string process_fastq_old(const string& fastq_path,
                 const vector<pair<string, string>>& adapters,
                 double threshold,
                 unsigned int thread_count,
                 const string& output_path) 
{
    
    // 线程共享资源
    queue<ReadData> read_queue;
    mutex queue_mutex;
    condition_variable cv;
    bool done_reading = false;
    size_t read_counter = 0;
    size_t segment_counter = 0;

    // 输出文件互斥锁
    mutex output_mutex;
    ofstream out_file(output_path);
    ofstream out_summary(output_path + ".summary.tsv");

    // 工作线程函数
    auto worker = [&]() {
        while(true) {
            unique_lock<mutex> lock(queue_mutex);
            cv.wait(lock, [&]{ return !read_queue.empty() || done_reading; });
            
            if(read_queue.empty() && done_reading) return;
            
            if(!read_queue.empty()) {
                ReadData data = read_queue.front();
                read_queue.pop();
                lock.unlock();

                // 执行比对和分段
                auto matches = find_adapters(data.sequence, adapters, threshold);
                if (matches.empty()) {
                    lock_guard<mutex> out_lock(output_mutex);
                    out_file << "@" << data.header << "|Noadapter" <<  "\n"
                        << data.sequence << "\n+\n" << data.quality << "\n";

                    out_summary << data.header << "\t" << "Noadapter" << "\t"
                        << "0" << "\t" << to_string(data.sequence.size()) << "\t" << "*" << "\n";
                    segment_counter++;
                }else {
                    auto segments = generate_segments(matches, data.header, data.sequence.length());
                    lock_guard<mutex> out_lock(output_mutex);
                    for(const auto& seg : segments) {
                        int start_idx = seg.start - 1;
                        int length = seg.end - seg.start + 1;
                        if(length < 100) continue;
                        if(start_idx < 0 || start_idx + length > data.sequence.size()) continue;
                        string subseq,subqual;
                        if(seg.strand == "+"){
                            int reslen = data.sequence.size()-start_idx-length;
                            int extendlen = std::min(20, reslen); // extend for anchor to extract CID
                            subseq = data.sequence.substr(start_idx, length+extendlen);
                            subqual = data.quality.substr(start_idx, length+extendlen);
                        }
                        
                        if(seg.strand == "-") {
                            int extendlen = std::min(20, start_idx);
                            subseq = data.sequence.substr(start_idx-extendlen, length+ extendlen);
                            subqual = data.quality.substr(start_idx-extendlen, length+ extendlen);
                            subseq = reverse_complement(subseq);
                            reverse(subqual.begin(), subqual.end());
                        }
                        
                        out_file << "@" << seg.read_id << "|" << seg.adapter_id << "|"
                                << seg.start << "-" << seg.end << "(" << seg.strand << ")\n"
                                << subseq << "\n+\n" << subqual << "\n";
                        
                        out_summary << seg.read_id << "\t" << seg.adapter_id << "\t"
                                  << seg.start << "\t" << seg.end << "\t" << seg.strand << "\n";
                        segment_counter++;
                    }
                }
            }
        }
    };

    // 启动工作线程
    vector<thread> workers;
    for(unsigned int t = 0; t < thread_count; ++t) {
        workers.emplace_back(worker);
    }

    // 主线程读取数据
    ifstream in_file(fastq_path);
    string line;
    while(getline(in_file, line)) {
        if(line.empty()) continue;
        
        ReadData data;
        data.header = "process_" + to_string(++read_counter);
        
        // 读取序列
        getline(in_file, data.sequence);
        while(in_file.peek() != '+' && in_file.peek() != EOF) {
            string append_seq;
            getline(in_file, append_seq);
            data.sequence += append_seq;
        }
        
        // 跳过+
        in_file.ignore(numeric_limits<streamsize>::max(), '\n');
        
        // 读取质量
        getline(in_file, data.quality);
        while(in_file.peek() != '@' && in_file.peek() != EOF) {
            string append_qual;
            getline(in_file, append_qual);
            data.quality += append_qual;
        }
        
        // 加入队列
        {
            lock_guard<mutex> lock(queue_mutex);
            read_queue.push(data);
        }
        cv.notify_one();
    }

    // 结束标志
    {
        lock_guard<mutex> lock(queue_mutex);
        done_reading = true;
    }
    cv.notify_all();

    // 等待线程结束
    for(auto& worker : workers) {
        worker.join();
    }

    string summary;
    summary = "Total reads: " + to_string(read_counter) + "\t" + "Processed segments: " + to_string(segment_counter) + "\n";
    return summary;
}



string fastq_segment_main(const string& fastq_path,
    const string& adapter_path,
    double threshold,
    unsigned int thread_count,
    const string& output_path) {
    auto adapters = load_adapters(adapter_path);
    // 添加锚点序列
    string anchor_prr = "CTGATAAGGTCGCCATGCCT";
    string anchor_flw = "TCTGCTGACGTACTGAGAGGCG";

    for(size_t i = 0; i < adapters.size(); ++i) {
        string& seq = adapters[i].second;
        if(i == 0) {
            seq += anchor_flw;
        } else if(i == adapters.size() - 1) {
            seq = anchor_prr + seq;
        } else {
            seq = anchor_prr + seq + anchor_flw;
        }
    }
    precompute_rc(adapters);
    string summary;
    summary = process_fastq(fastq_path, adapters, threshold, thread_count, output_path);
    return summary;
}
} // namespace FastqProcessor

