#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <queue>
#include <condition_variable>
#include <htslib/sam.h>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <set>
#include "process_stereo_seq_bam.h"

using namespace std;

namespace HandleStereoBam {
// 线程安全输出器
class ConcurrentWriter {
private:
    mutex mtx;
    ofstream out;
    atomic<size_t> counter{0};
    static const size_t BATCH_SIZE = 1000;
    vector<string> buffer;

public:
    ConcurrentWriter(const string& filename) : out(filename) {}

    void write(const string& line) {
        lock_guard<mutex> lock(mtx);
        buffer.push_back(line);
        if (buffer.size() >= BATCH_SIZE) {
            flush();
        }
    }

    void flush() {
        for (const auto& line : buffer) {
            out << line << "\n";
            ++counter;
        }
        buffer.clear();
        out.flush();
    }

    ~ConcurrentWriter() {
        lock_guard<mutex> lock(mtx);
        flush();
    }

    size_t count() const { return counter.load(); }
};

struct GeneRange {
    string name;
    string region;
};
// 高效分割函数
vector<string> split(const string& s, const string& delim) {
    vector<string> tokens;
    size_t pos = 0, last = 0;
    while ((pos = s.find(delim, last)) != string::npos) {
        if (pos > last)
            tokens.push_back(s.substr(last, pos - last));
        last = pos + delim.length();
    }
    if (last < s.length())
        tokens.push_back(s.substr(last));
    return tokens;
}

void save_cidpgene(const unordered_map<string, vector<string>>& cidpgene, const string& filename) {
    ofstream out(filename);
    for (const auto& entry : cidpgene) {
        out << entry.first;
        for (const auto& gene : entry.second) {
            out << "\t" << gene;
        }
        out << "\n";
    }
    out.close();
    cout << "Saved cidpgene to " << filename << " with " << cidpgene.size() << " entries" << endl;
}

void load_cidpgene(unordered_map<string, vector<string>>& cidpgene, const string& filename) {
    ifstream in(filename);
    string line;
    while (getline(in, line)) {
        istringstream ss(line);
        string pos;
        ss >> pos;
        vector<string> genes;
        string gene;
        while (ss >> gene) {
            genes.push_back(gene);
        }
        // 去重
        set<string> unique_genes(genes.begin(), genes.end());
        cidpgene[pos] = vector<string>(unique_genes.begin(), unique_genes.end());
    }
    in.close();
    cout << "Loaded cidpgene from " << filename << " with " << cidpgene.size() << " entries" << endl;
}

vector<GeneRange> load_gene_ranges(const string& gtf_file) {
    vector<GeneRange> ranges;
    ifstream file(gtf_file);
    string line;

    while (getline(file, line)) {
        if (line[0] == '#') continue;

        istringstream ss(line);
        string seqname, feature, start_str, end_str, attributes;
        ss >> seqname >> feature >> feature; // 跳过source和feature

        if (feature != "gene") continue;

        ss >> start_str >> end_str;
        uint64_t start = stoull(start_str);
        uint64_t end = stoull(end_str);

        // 扩展区域
        start = start > 1000 ? start - 1000 : 0;
        end += 1000;

        // 解析gene_name
        string gene_name;
        getline(ss, attributes);
        size_t pos = attributes.find("gene_id");
        if (pos != string::npos) {
            size_t qstart = attributes.find('"', pos) + 1;
            size_t qend = attributes.find('"', qstart);
            gene_name = attributes.substr(qstart, qend - qstart);
        }

        if (!gene_name.empty()) {
            ranges.push_back({gene_name, seqname + ":" + to_string(start) + "-" + to_string(end)});
        }
    }
    return ranges;
}

void bam_cidpos(const GeneRange& gene,
                 const string& bam_path,
                unordered_map<string, vector<string>>& cidpgene,
            mutex &mtx) {
    samFile* bam = sam_open(bam_path.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(bam);
    hts_idx_t* idx = sam_index_load(bam, bam_path.c_str());
    bam1_t* rec = bam_init1();

    hts_itr_t* iter = sam_itr_querys(idx, header, gene.region.c_str());

    while (sam_itr_next(bam, iter, rec) >= 0) {
        const string qname = bam_get_qname(rec);
        // 提取pos信息
        auto segments = split(qname, "|||");
        string pos = "";
        if (segments.size() < 2) continue;
        for(auto& segment : segments){
            if(segment.substr(0,5) == "CB:Z:"){
                pos = segment.substr(5);
            }
        }
        mtx.lock();
        cidpgene[pos].push_back(gene.name);
        mtx.unlock();
    }
    hts_itr_destroy(iter);
    bam_destroy1(rec);
    hts_idx_destroy(idx);
    sam_close(bam);
    bam_hdr_destroy(header);
}

void bam_cidpos_tag(const GeneRange& gene,
                 const string& bam_path,
                unordered_map<string, vector<string>>& cidpgene,
            mutex &mtx)
{
    samFile* bam = sam_open(bam_path.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(bam);
    hts_idx_t* idx = sam_index_load(bam, bam_path.c_str());
    if (!idx) {  // 检查索引是否加载成功
        std::cerr << "Failed to load BAM index for " << bam_path << endl;
    }

    bam1_t* rec = bam_init1();
    //cout << gene.region.c_str() << endl;
    hts_itr_t* iter = sam_itr_querys(idx, header, gene.region.c_str());
    int cx = -1, cy = -1;

    int count = 0;
    while (sam_itr_next(bam, iter, rec) >= 0) {
        //const string qname = bam_get_qname(rec);
        //if (rec->core.flag & BAM_FUNMAP) continue; // 跳过未映射记录
        int cx = -1, cy = -1;
        uint8_t *cx_ptr = bam_aux_get(rec, "Cx");
        uint8_t *cy_ptr = bam_aux_get(rec, "Cy");

        if (cx_ptr && cy_ptr) {
            cx = bam_aux2i(cx_ptr);
            cy = bam_aux2i(cy_ptr);
            string pos = to_string(cx) + "_" + to_string(cy);
            mtx.lock();
            cidpgene[pos].push_back(gene.name);
            mtx.unlock();
            count++;
        }
    }
    cout << "gene record:" << count << endl;
    //cout << "cidpos" << cidpgene.size() << endl;
    hts_itr_destroy(iter);
    bam_destroy1(rec);
    hts_idx_destroy(idx);
    sam_close(bam);
    bam_hdr_destroy(header);
}


void process_range_chunk(const vector<GeneRange>& genes,
                        string tag,
                        size_t start,
                        size_t end,
                        const string& bam_path,
                        unordered_map<string, vector<string>>& cidpgene,
                    mutex &mtx) {
    if(tag == "T"){
        for (size_t i = start; i <= end; ++i) {
            bam_cidpos_tag(genes[i], bam_path, cidpgene,mtx);
        }
    }else{
        for (size_t i = start; i <= end; ++i) {
            bam_cidpos(genes[i], bam_path, cidpgene,mtx);
        }
    }

}

void process_bam_extract(string gtf_file, string bam_file, int num_threads,string tag,
                         unordered_map <string, vector<string>> &cidpgene) {
    // 加载基因区域
    auto gene_ranges = load_gene_ranges(gtf_file);
    if (gene_ranges.empty()) {
        std::cerr << "No valid gene ranges found" << endl;
    }

    // 分配任务
    const size_t total_genes = gene_ranges.size();
    const size_t chunk_size = (total_genes + num_threads - 1) / num_threads;
    mutex mtx;
    vector<thread> workers;
    for (int i = 0; i < num_threads; ++i) {
        const size_t start = i * chunk_size;
        const size_t end = min(start + chunk_size - 1, total_genes - 1);

        if (start >= total_genes) break;
        workers.emplace_back([&, start, end]() {
            process_range_chunk(gene_ranges, tag, start, end, bam_file, ref(cidpgene),ref(mtx));
        });
    }
    // 等待所有线程完成
    for (auto& t : workers) {
        t.join();
    }
}


// 构建参考数据哈希表
void build_reflist(const string& filename, unordered_map<string, string>& refmap) {
    ifstream file(filename);
    string line;

    while (getline(file, line)) {
        auto parts = split(line, "\t");
        if (parts.size() < 3) continue;

        string pos = parts[1] + "_" + parts[2];
        refmap[pos] = parts[0] + "\t" + parts[1] + "\t" + parts[2];
    }
    cout << "Loaded " << refmap.size() << " reference positions" << endl;
}

// 处理基因数据文件
void Joint_cidseq(unordered_map <string, vector<string>> & cidpgene, unordered_map<string, string>& refmap,
                    const string& outfile) {
    ofstream output(outfile);
    string line;
    size_t count = 0;

    for(auto cidg : cidpgene){
        string pos = cidg.first;           // 键（string 类型）
        vector<string> genes = cidg.second;
        unordered_set<string> unique_set(genes.begin(), genes.end());
        genes.assign(unique_set.begin(), unique_set.end());

        auto it = refmap.find(pos);
        if (it != refmap.end()) {
            for (const auto& gene : genes) {
                output << refmap[pos] << "\t" << gene << "\n";
                if (++count % 1000000 == 0)
                    cout << "Processed " << count << " records" << endl;
            }
        }
    }
    cout << "Total merged records: " << count << endl;
}


void merge_cidlist_old(string refcid_file, unordered_map <string, vector<string>> &cidpgene, string out_file) {
    // 内存优化：设置更大哈希表桶数量
    unordered_map<string, string> cidseq;
    //cidseq.reserve(150000000);  // 预分配1.5亿桶
    // 第一步：加载参考数据
    build_reflist(refcid_file, cidseq);
    // 第二步：处理基因数据
    Joint_cidseq(cidpgene, cidseq, out_file);
}

void build_reftable_old(string gtf_file, string bam_file, string refcid_file, string outfile, int num_threads,string tag) {

    unordered_map <string, vector<string>> cidpgene;
    cout << "extract from bam..." << endl;
    process_bam_extract(gtf_file, bam_file, num_threads, tag, cidpgene);
    cout << "cid number:" << cidpgene.size() << endl;
    string outfile2 = outfile + ".gene_coord_CID.tsv";
    cout << "Merge with cid seq..." << endl;
    merge_cidlist_old(refcid_file, cidpgene, outfile2);
}

void merge_cidlist(string refcid_file, unordered_map<string, vector<string>>& cidpgene, string out_file) {
    // 第一步：先保存cidpgene到临时文件
    string cidpgene_file = out_file + ".cidpgene.tmp";
    save_cidpgene(cidpgene, cidpgene_file);

    // 清空原始cidpgene释放内存
    cidpgene.clear();

    // 第二步：分批处理ref文件，避免内存不足
    ifstream ref_file(refcid_file);
    ofstream output(out_file);
    string line;
    size_t batch_size = 1000000; // 每批处理100万行
    size_t processed_count = 0;
    size_t output_count = 0;

    cout << "Starting merge process with batch size: " << batch_size << endl;

    // 预分配内存优化[8,9](@ref)
    unordered_map<string, vector<string>> batch_cidpgene;
    batch_cidpgene.reserve(min(batch_size, cidpgene.size()));

    // 重新加载cidpgene（如果需要可以分批加载）
    load_cidpgene(cidpgene, cidpgene_file);

    // 设置合适的负载因子以减少哈希冲突[9](@ref)
    cidpgene.max_load_factor(0.75);

    // 分批读取ref文件并处理
    vector<string> ref_batch;
    ref_batch.reserve(batch_size);

    while (getline(ref_file, line)) {
        ref_batch.push_back(line);
        processed_count++;

        if (ref_batch.size() >= batch_size) {
            // 处理当前批次
            for (const auto& ref_line : ref_batch) {
                auto parts = split(ref_line, "\t");
                if (parts.size() < 3) continue;

                string pos = parts[1] + "_" + parts[2];

                auto it = cidpgene.find(pos);
                if (it != cidpgene.end()) {
                    for (const auto& gene : it->second) {
                        output << parts[0] << "\t" << parts[1] << "\t" << parts[2] << "\t" << gene << "\n";
                        output_count++;
                    }
                }
            }

            ref_batch.clear();
            ref_batch.reserve(batch_size);

            if (processed_count % (batch_size * 10) == 0) {
                cout << "Processed " << processed_count << " ref records, output " << output_count << " records" << endl;
            }
        }
    }

    // 处理最后一批
    if (!ref_batch.empty()) {
        for (const auto& ref_line : ref_batch) {
            auto parts = split(ref_line, "\t");
            if (parts.size() < 3) continue;

            string pos = parts[1] + "_" + parts[2];

            auto it = cidpgene.find(pos);
            if (it != cidpgene.end()) {
                for (const auto& gene : it->second) {
                    output << parts[0] << "\t" << parts[1] << "\t" << parts[2] << "\t" << gene << "\n";
                    output_count++;
                }
            }
        }
    }

    ref_file.close();
    output.close();

    // 清理临时文件
    remove(cidpgene_file.c_str());

    cout << "Merge completed. Total processed: " << processed_count << ", output: " << output_count << " records" << endl;
}

// 进一步优化的版本：使用流式处理，完全避免大内存占用
void merge_cidlist_streaming(string refcid_file, unordered_map<string, vector<string>>& cidpgene, string out_file) {
    // 保存cidpgene到文件
    string cidpgene_file = out_file + ".cidpgene.tmp";
    save_cidpgene(cidpgene, cidpgene_file);

    // 使用更内存高效的方式：逐行处理ref文件
    ifstream ref_file(refcid_file);
    ofstream output(out_file);
    string line;
    size_t processed_count = 0;
    size_t output_count = 0;

    cout << "Starting streaming merge process..." << endl;

    // 逐行处理，避免批量存储
    while (getline(ref_file, line)) {
        auto parts = split(line, "\t");
        if (parts.size() < 3) continue;

        string pos = parts[1] + "_" + parts[2];

        // 直接在cidpgene中查找
        auto it = cidpgene.find(pos);
        if (it != cidpgene.end()) {
            for (const auto& gene : it->second) {
                output << parts[0] << "\t" << parts[1] << "\t" << parts[2] << "\t" << gene << "\n";
                output_count++;
            }
        }

        processed_count++;

        if (processed_count % 1000000 == 0) {
            cout << "Processed " << processed_count << " ref records, output " << output_count << " records" << endl;
            // 定期刷新输出缓冲区[6](@ref)
            output.flush();
        }
    }

    ref_file.close();
    output.close();

    // 清理临时文件
    remove(cidpgene_file.c_str());

    cout << "Streaming merge completed. Total processed: " << processed_count << ", output: " << output_count << " records" << endl;
}


void build_reftable(string gtf_file, string bam_file, string refcid_file, string outfile, int num_threads, string tag) {
    unordered_map<string, vector<string>> cidpgene;

    cout << "Extracting from bam..." << endl;
    process_bam_extract(gtf_file, bam_file, num_threads, tag, cidpgene);
    cout << "CID number: " << cidpgene.size() << endl;
    string outfile2 = outfile + ".gene_coord_CID.tsv";
    cout << "Merge with CID seq..." << endl;

    // 根据cidpgene的大小选择合适的合并策略
    if (cidpgene.size() > 1000000) {
        // 对于大数据集，使用流式处理
        cout << "Using streaming merge for large dataset..." << endl;
        merge_cidlist_streaming(refcid_file, cidpgene, outfile2);
    } else {
        // 对于小数据集，使用批量处理
        cout << "Using batch merge for small dataset..." << endl;
        merge_cidlist(refcid_file, cidpgene, outfile2);
    }
}

}
