# Fullscope-seq

Project Description

A C++ toolset for processing Long-reads Spatial Transcriptomic Sequencing data (Stereo-Seq) , providing a complete analysis pipeline from raw FASTQ files to CID (Cell Identifier) mapping.


# System Requirements

• Linux operating system

• C++17 compatible compiler

• minimap2, samtools (installed and available in PATH)

• Sufficient disk space for intermediate files

# Installation

## Clone the repository
`git clone (https://github.com/wwei-lab/Fullscope`

# Quick Start

Complete Analysis Pipeline

## Run complete analysis pipeline (from FASTQ to CID mapping)
`./Fullscope-1.5 count input.fq adapters.fa anchor.fa 0.8 genes.gtf genome.fa index.cid index_threshold.txt 7 6 16 ./output sample1`

# Step-by-Step Analysis

## 1. FASTQ Sequence Segmentation

### Segment sequences in FASTQ file
`./Fullscope-1.5 process_fq input.fq adapters.fa anchor.fa 0.8 16 segmented.fq`


## 2. Sequence Alignment (minimap2)

### Align using minimap2 (requires separate installation)
`minimap2 -K500m --secondary=no -a -x splice --splice-flank=yes -t 16 genome.fa segmented.fq | samtools sort > aligned.bam`
`samtools index aligned.bam`


## 3. CID Extraction

### Extract CID from BAM file
`./Fullscope-1.5 extract aligned.bam genes.gtf extracted_cid.tsv 16`


## 4. CID Mapping

### Map CID to index
`./Fullscope-1.5 map extracted_cid.tsv index.cid index_threshold.txt mapped_cid.tsv 16 7 6`


# Detailed Parameter Description

## count Command (Complete Pipeline)

`./Fullscope-1.5 count <input.fq> <adapters.fa> <anchor.fa> <segment_threshold> <refgtf> <refgenome> <index.cid> <index_threshold.txt> <kmer> <bucketnum> <threads> <outfold> <outprex>`


## Parameters:
• input.fq: Input FASTQ file

• adapters.fa: Adapter sequence file

• anchor.fa: Anchor sequence file

• segment_threshold: Segmentation threshold (0-1)

• refgtf: Reference genome GTF file

• refgenome: Reference genome FASTA file

• index.cid: CID index file

• index_threshold.txt: Index threshold file

• kmer: k-mer length (default: 7)

• bucketnum: Number of buckets (default: 6)

• threads: Number of threads

• outfold: Output folder

• outprex: Output file prefix

## Other Available Commands

process_fq - FASTQ Processing

`./Fullscope-1.5 process_fq input.fq adapters.fa anchor.fa 0.8 16 output.fq`


extract - CID Extraction from BAM

`./Fullscope-1.5 extract input.bam genes.gtf output.tsv 16`


extract_fq - CID Extraction from FASTQ

`./Fullscope-1.5 extract_fq input.fq output.tsv 16`


map - CID Mapping

`./Fullscope-1.5 map reads.txt index.cid index_threshold.txt output.txt 16 7 6`


map_p - Precise CID Mapping

`./Fullscope-1.5 map_p reads.txt index.cid output.txt 16 7`

build_idx - Index Building

 Build fast index
`./Fullscope-1.5 build_idx f input.txt 16 7 6 index.cid`

 Build precise index
`./Fullscope-1.5 build_idx p input.txt 16 7 index_precise.cid`


bamtoref - BAM to Reference Table

`./Fullscope-1.5 bamtoref genes.gtf input.bam ref_list.txt output 16 T`


# Output Structure

Complete pipeline generates the following directory structure:

output/
├── Segment/           # Segmented FASTQ files
├── Alignment/         # Alignment BAM files
├── CIDextract/        # Extracted CID files
├── CIDmap/           # Mapped CID files
└── sample1.summary.txt # Analysis summary


# Dependencies

• minimap2: Required for sequence alignment

• samtools: Required for BAM file processing

• C++17: Compiler with C++17 support

# Notes

1. Dependencies: Ensure minimap2 and samtools are properly installed and accessible in PATH
2. Memory Requirements: Large datasets may require significant memory (recommended ≥32GB)
3. Thread Configuration: Set thread count appropriately based on server configuration
4. File Permissions: Ensure sufficient disk space and write permissions

# Troubleshooting

Common Issues

1. Compilation Errors: Check if g++ version supports C++17
2. Runtime Errors: Verify all input files exist and have correct formats
3. Insufficient Memory: Reduce thread count or increase system memory

# Getting Help

For help, use:
./Fullscope-1.5


# Version Information

• Current Version: 1.0

• Last Updated: January 10, 2026

• Supported Systems: Linux

# License

Apache License 2.0

This project is licensed under the Apache License, Version 2.0. See the LICENSE file for details.
