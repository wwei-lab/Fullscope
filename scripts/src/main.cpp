#include "CID_extract.h"
#include "CID_mapping.h"
#include "fastq_segment.h"
#include "process_stereo_seq_bam.h"
#include "CID_index.h"
#include "CID_index_precise.h"
#include "CID_mapping_precise.h"
#include "CID_extract_nogene.h"
#include "utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <filesystem>
using namespace std;


void print_usage() {
    cerr << "Usage:\n"
         << "  count <input.fq> <adapters.fa> <segment_threshold> <refgtf> <refgenome> <index.cid> <index_threshold.txt> <kmer=7> <bucketnum=6> <threads> <outfold> <outprex>\n"
         << "  process_fq <input.fq> <adapters.fa> <threshold> <threads> <output.fq>\n"
         << "  extract   <input.bam> <genes.gtf> <output.tsv> <threads=4> \n"
         << "  extract_fq   <input.fq> <output.tsv> <threads=4> \n"
         << "  map       <reads.txt> <index.cid> <index_threshold.txt> <output.txt> <threads=4> <kmer=7> <bucketnum=6>\n"
         << "  map_p       <reads.txt> <index.cid> <output.txt> <threads=4> <kmer=7>\n"
         << "  build_idx <type> <input.txt> <threads> <kmer> <bucketnum> <output>\n"
         << "  bamtoref <gtf> <bam> <reflist> <output> <threads=4> <tag=T>\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage();
        return 1;
    }

    const string command = argv[1];

    try {
        if (command == "count") {
            if (argc < 13) {
                cerr << "extract requires at least 11 arguments\n";
                print_usage();
                return 1;
            }
            string input_fq = argv[2];
            string adapter_fa = argv[3];
            double seg_thred = stod(argv[4]);
            string gtf = argv[5];
            string genome = argv[6];
            string index = argv[7];
            string index_thred = argv[8];
            int kmer = stoi(argv[9]);
            int bucketnum = stoi(argv[10]);
            int threads = stoi(argv[11]);
            string outfold = argv[12];
            string outprex = argv[13];
            
            string summary1;
            // fastq segment #################
            cout << "processing segmentation..." << endl;
            string outfold1 = outfold + "/Segment/";
            if (!std::filesystem::exists(outfold1)) {
                std::filesystem::create_directories(outfold1);
            } else {
                std::cerr << "failed to create directory" << std::endl;
                return -1;
            }
            string outfile1 = outfold1 + outprex + ".fq";
            summary1 = FsFastqSegment::fastq_segment_main(input_fq, adapter_fa, seg_thred, threads, outfile1);

            // minimap2 ##################
            cout << "processing minimap2 mapping..." << endl;
            string outfold2 = outfold + "/Alignment/";
            if (!std::filesystem::exists(outfold2)) {
                std::filesystem::create_directories(outfold2);
            }
            string outfile2 = outfold2 + outprex + ".bam";
            string command1 = "minimap2 -K500m --secondary=no -a -x splice --splice-flank=yes -t " + to_string(threads) + " " + 
            genome + " " + outfile1 + " | samtools sort " + "> " + outfile2;
            int ret = system(command1.c_str());
            if (ret != 0) {
                std::cerr << "failed to process minimap2" << std::endl;
                return -1;
            }
            string command2 = "samtools index " + outfile2;
            ret = system(command2.c_str());
            if (ret != 0) {
                std::cerr << "failed to build bam index" << std::endl;
                return -1;
            }

            // CID extracting ##################
            cout << "processing CID extracting..." << endl;
            string outfold3 = outfold + "/CIDextract/";
            string outfile3 = outfold3 + outprex + "_raw_CID.tsv";
            if (!std::filesystem::exists(outfold3)) {
                std::filesystem::create_directories(outfold3);
            }
            CIDExtract::GeneRanges gene_ranges;
            CIDExtract::load_gene_ranges(gtf, gene_ranges);
            CIDExtract::Config config;
            config.num_threads = threads;
            string summary2;
            summary2 = CIDExtract::process_bam(outfile2, outfile3, gene_ranges, config);

            // CID mapping ##################
            cout << "processing CID mapping..." << endl;
            string outfold4 = outfold + "/CIDmap/";
            string outfile4 = outfold4 + outprex + "_mapped_CID.tsv";
            if (!std::filesystem::exists(outfold4)) {
                std::filesystem::create_directories(outfold4);
            }
            string summary3;
            summary3 = CIDmapping::CIDmap_processer(outfile3, index, threads, kmer, bucketnum, outfile4, index_thred);
            cout << "All processes were successfully done!" << endl;
            string outfiles = outfold + outprex + ".summary.txt";
            ofstream out(outfiles);
            out << summary1 << summary2 << summary3;
            
        } else if (command == "process_fq") {
            if (argc < 7) {
                cerr << "process_fq requires 6 arguments\n";
                print_usage();
                return 1;
            }
            double threshold = stod(argv[4]);
            int threads = stoi(argv[5]);
            string summary = FsFastqSegment::fastq_segment_main(argv[2], argv[3], threshold, threads, argv[6]);
            cout << summary << endl;

        } else if (command == "extract") {
            if (argc < 5) {
                cerr << "extract requires at least 4 arguments\n";
                print_usage();
                return 1;
            }
            
            CIDExtract::GeneRanges gene_ranges;
            CIDExtract::load_gene_ranges(argv[3], gene_ranges);
            
            CIDExtract::Config config;
            config.num_threads = stoi(argv[5]);
            
            CIDExtract::process_bam(argv[2], argv[4], gene_ranges, config);

        } else if (command == "extract_fq") {
            if (argc < 4) {
                cerr << "extract requires at least 3 arguments\n";
                print_usage();
                return 1;
            }
            string inputfq = argv[2];
            string output = argv[3];
            CIDExtractFastq::Config config;
            config.num_threads = stoi(argv[4]);
            
            CIDExtractFastq::extractcid_fastq(inputfq, output, config);

        } else if (command == "map") {
            if (argc < 7) {
                cerr << "map requires at least 5 arguments\n";
                print_usage();
                return 1;
            }
            int threads,kmer,bucketnum;

            if (argc >= 7) threads = stoi(argv[6]);
            if (argc >= 8) kmer = stoi(argv[7]);
            if (argc >= 9) bucketnum = stoi(argv[8]);

            CIDmapping::CIDmap_processer(argv[2], argv[3], threads, kmer, bucketnum,argv[5],argv[4]);

        } else if (command == "map_p") {
            if (argc < 6) {
                cerr << "map requires at least 4 arguments\n";
                print_usage();
                return 1;
            }
            int threads,kmer,bucketnum;

            threads = stoi(argv[5]);
            kmer = stoi(argv[6]);
            string readfile = argv[2];
            string indexfile = argv[3];
            string outfile = argv[4];
            CIDmappingPrecise::CIDmap_processer_precise(readfile, indexfile, threads, kmer, outfile);

        } else if (command == "build_idx") {
            if (argc < 7) {
                cerr << "build_idx requires 5 arguments\n";
                print_usage();
                return 1;
            }
            string type = argv[2];
            string cidref = argv[3];
            int threads = stoi(argv[4]);
            int kmerlen = stoi(argv[5]);
            int bucketnum = stoi(argv[6]);
            string outfile = argv[7];

            if(type == "f"){
                CIDindex::build_index(cidref, kmerlen, bucketnum, threads, outfile);
            }else if(type == "p"){
                CIDindexPrecise::build_index_precise(cidref, kmerlen, threads, outfile);
            }

        } else if (command == "bamtoref") {
            // "  bamtoref <gtf> <bam> <reflist> <output> <threads=4>\n"
            if (argc < 8) {
                cerr << "bamtoref requires 6 arguments\n";
                print_usage();
                return 1;
            }
            string gtf_file = argv[2];
            string ref_bam = argv[3];
            string ref_cid = argv[4];
            string outprex = argv[5];
            int threads = stoi(argv[6]);
            string tag = argv[7];
            //string gtf_file, string bam_file, string refcid, string outfile, int num_threads
            HandleStereoBam::build_reftable(gtf_file, ref_bam, ref_cid, outprex, threads, tag);

        } else {
            cerr << "Unknown command: " << command << "\n";
            return 1;
        }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}