# Fullscope-seq  
**A C++ toolset for  long-read spatial transcriptomic data (Stereo-seq with ONT/Pacbio/Cyclone): from raw FASTQ to CID mapping in one command.**

---

##  Overview  
Fullscope-seq is a high-performance pipeline that turns **Stereo-seq long-read spatial transcriptomic FASTQ** into **cell-resolution CID maps**.  
Everything is written in C++17 and engineered for Linux servers with ≥32 GB RAM.

---

##  Quick Start (5 commands)

1. Clone  
   ```bash
   git clone https://github.com/wwei-lab/Fullscope.git
   cd Fullscope
   ```

2. Run the complete pipeline  
   ```bash
   ./Fullscope-1.5 count \
     input.fq adapters.fa anchor.fa 0.8 \
     genes.gtf genome.fa index.cid index_threshold.txt \
     7 6 16 ./output sample1
   ```

3. Done – results are in `output/`  
   ```
   output/
   ├── Segment/           # demultiplexed reads
   ├── Alignment/         # sorted & indexed BAM
   ├── CIDextract/        # per-read CID tables
   ├── CIDmap/            # CID → spatial index
   └── sample1.summary.txt
   ```

---

## System Requirements

| Component       | Requirement               |
|-----------------|---------------------------|
| OS              | Linux (kernel ≥3.10)      |
| Compiler        | g++ ≥7 / clang++ ≥5 (C++17) |
| Third-party     | minimap2, samtools in `$PATH` |
| Hardware        | ≥32 GB RAM, ≥16 cores recommended |

---

## Installation

```bash
# 1. clone
git clone https://github.com/wwei-lab/Fullscope.git && cd Fullscope

# 2. build (CMake ≥3.14)
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
# binary is now ./Fullscope-1.5
```

---

## Step-by-step (if you prefer modular runs)

| Step | Command |
|------|---------|
| 1. FASTQ segmentation | `./Fullscope-1.5 process_fq input.fq adapters.fa anchor.fa 0.8 16 segmented.fq` |
| 2. Alignment | `minimap2 -K500m --secondary=no -a -x splice --splice-flank=yes -t 16 genome.fa segmented.fq \| samtools sort -o aligned.bam` |
| 3. CID extraction | `./Fullscope-1.5 extract aligned.bam genes.gtf extracted_cid.tsv 16` |
| 4. CID mapping | `./Fullscope-1.5 map extracted_cid.tsv index.cid index_threshold.txt mapped_cid.tsv 16 7 6` |

---

## Parameter Bible

### `count` – single-command pipeline
```
./Fullscope-1.5 count \
  <input.fq> <adapters.fa> <anchor.fa> <seg_thresh> \
  <ref.gtf> <ref.fa> <index.cid> <index_thresh.txt> \
  <k> <buckets> <threads> <outDir> <prefix>
```

| Param | Meaning | Example |
|-------|---------|---------|
| `seg_thresh` | segment confidence (0–1) | `0.8` |
| `k` | k-mer length | `7` |
| `buckets` | hash buckets | `6` |
| `threads` | CPU cores | `16` |

---

### Auxiliary tools

| Tool | One-liner |
|------|-----------|
| `process_fq` | `./Fullscope-1.5 process_fq in.fq adapters.fa anchor.fa 0.8 16 out.fq` |
| `extract` | `./Fullscope-1.5 extract aligned.bam genes.gtf cid.tsv 16` |
| `extract_fq` | `./Fullscope-1.5 extract_fq in.fq cid.tsv 16` |
| `map` | `./Fullscope-1.5 map cid.tsv index.cid thresh.txt out.tsv 16 7 6` |
| `map_p` (precise) | `./Fullscope-1.5 map_p cid.tsv index.cid out.tsv 16 7` |
| `build_idx` | fast: `./Fullscope-1.5 build_idx f list.txt 16 7 6 index.cid`  
precise: `./Fullscope-1.5 build_idx p list.txt 16 7 index_precise.cid` |
| `bamtoref` | `./Fullscope-1.5 bamtoref genes.gtf in.bam ref_list.txt out 16 T` |

---

## Output Structure

```
output/
├── Segment/           # demultiplexed *.fq
├── Alignment/         # *.bam + *.bai
├── CIDextract/        # *.cid.tsv
├── CIDmap/            # *.mapped.cid.tsv
└── <prefix>.summary.txt
```

---

## Notes & Pro-tips

1. **Dependencies** – ensure `minimap2` and `samtools` are in `$PATH`.  
2. **Memory** – 1 GB per 1 M reads is a safe rule of thumb.  
3. **Threads** – leave 2 cores free for I/O.  
4. **Disk** – reserve 5× the input FASTQ size for intermediates.

---

## Help & Version

```bash
./Fullscope-1.5          # print help
./Fullscope-1.5 --version # 1.5.0 (Jan 2026)
```

