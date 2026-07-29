# Fullscope

Fullscope is the reference software for **Fullscope-seq**, a full-length,
single-molecule, large-field-of-view spatial transcriptomics method. The
software processes Stereo-seq-compatible long-read data from concatenated
cDNA segmentation through spatial CID assignment and optional transcript
annotation.

Associated paper:

> Liu H, Hong Y, Zhang YS, *et al.* Full-length single-cell spatial
> transcriptomics reveals spatial and cell-type-specific transcript isoforms
> in the primate brain. *Nature Methods* (2026).
> [https://doi.org/10.1038/s41592-026-03174-y](https://doi.org/10.1038/s41592-026-03174-y)

## What is included

- A C++23 core for:
  - programmed-concatemer FASTQ segmentation;
  - CID extraction from FASTQ or BAM;
  - fast and precise CID-index construction;
  - error-tolerant CID mapping.
- Portable command-line wrappers for:
  - FASTQ/BAM input;
  - Stereo-seq barcode-index preparation;
  - splice-aware alignment;
  - optional Bambu transcript annotation;
  - merging transcript assignments with spatial coordinates;
  - single-sample and Slurm batch execution.
- Small test data and a smoke test.
- Analysis notebooks used for the associated study.

## Workflow

```text
concatenated long-read FASTQ or BAM
                  |
                  v
       cDNA segmentation (C++)
                  |
       +----------+-----------+
       |                      |
       v                      v
CID extraction/mapping   splice-aware alignment
       |                      |
       |                 optional Bambu
       +----------+-----------+
                  |
                  v
 spatially resolved transcript assignments
```

## Installation

### Recommended: conda or mamba

Fullscope is supported on 64-bit Linux. The build requires a C++23 compiler,
CMake, SeqAn3, cereal and HTSlib. The supplied core environment installs the
build and command-line dependencies; a complete environment with the optional
R/Bambu stack is also provided.

```bash
git clone https://github.com/wwei-lab/Fullscope.git
cd Fullscope

mamba env create -f environment-core.yml
conda activate fullscope

bash install.sh --prefix "$CONDA_PREFIX"
```

If `mamba` is unavailable, replace the first command with:

```bash
conda env create -f environment-core.yml
```

The core environment supports segmentation, CID processing, and alignment. To
also install the optional R/Bambu transcript-annotation dependencies, create
the complete environment instead:

```bash
mamba env create -f environment.yml
conda activate fullscope
bash install.sh --prefix "$CONDA_PREFIX"
```

Confirm the installation:

```bash
fullscope --help
fullscope segment --help
fullscope --version
```

### Build with an existing environment

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 8
cmake --install build --prefix "$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
```

## Smoke test

The bundled smoke test runs segmentation on 20 reads and does not require a
reference genome or Stereo-seq mask:

```bash
bash tests/smoke_test.sh "$(command -v fullscope)"
```

A successful run ends with `Smoke test passed`.

## Quick start

### 1. Segmentation only

Use this mode to validate installation or to split concatenated cDNA reads
before downstream processing:

```bash
fullscope segment \
  --raw-fq reads.fastq \
  --out results/sample_fragment.fastq \
  --threads 8
```

The packaged adapter and anchor FASTA files are used automatically. Override
them only when the library design differs:

```bash
fullscope segment \
  --raw-fq reads.fastq \
  --out results/sample_fragment.fastq \
  --adapter-fa custom_adapters.fa \
  --anchor-fa custom_anchor.fa \
  --segthreshold 0.15 \
  --threads 8
```

### 2. Complete workflow from FASTQ

Full CID-index creation requires
[ST_BarcodeMap](https://github.com/STOmics/ST_BarcodeMap), which is maintained
separately and is not bundled with Fullscope.

```bash
fullscope run \
  --sample sample01 \
  --outdir results/sample01 \
  --raw-fq reads.fastq.gz \
  --stereoindex sample01.barcodeToPos.h5 \
  --barcode-map /path/to/ST_BarcodeMap \
  --genome reference/genome.fa \
  --gtf reference/genes.gtf \
  --threads 32
```

Both uncompressed FASTQ and `.fastq.gz` input are accepted. Gzipped input is
decompressed into the sample output directory before C++ processing.

### 3. Complete workflow from BAM

```bash
fullscope run \
  --sample sample01 \
  --outdir results/sample01 \
  --input-bam reads.bam \
  --stereoindex sample01.barcodeToPos.h5 \
  --barcode-map /path/to/ST_BarcodeMap \
  --genome reference/genome.fa \
  --gtf reference/genes.gtf \
  --threads 32
```

### 4. Add transcript annotation and spatial merge

```bash
fullscope run \
  --sample sample01 \
  --outdir results/sample01 \
  --raw-fq reads.fastq.gz \
  --stereoindex sample01.barcodeToPos.h5 \
  --barcode-map /path/to/ST_BarcodeMap \
  --genome reference/genome.fa \
  --gtf reference/genes.gtf \
  --run-bambu \
  --merge-annot \
  --threads 32
```

## Configuration

For repeated runs, copy the portable example and fill in local reference paths:

```bash
cp fullscope_toolkit/config/config.example.env site.env
```

Then run:

```bash
fullscope run \
  --config-env site.env \
  --sample sample01 \
  --outdir results/sample01 \
  --raw-fq reads.fastq.gz \
  --stereoindex sample01.barcodeToPos.h5
```

Command-line flags take precedence over values loaded from the configuration
file.

## Main workflow options

| Option | Description |
|---|---|
| `--sample NAME` | Sample identifier; required. |
| `--outdir PATH` | Sample output directory; required. |
| `--raw-fq PATH` | Input FASTQ or FASTQ.GZ. |
| `--input-bam PATH` | Input BAM; used when FASTQ is not supplied. |
| `--stereoindex PATH` | Stereo-seq `barcodeToPos.h5` mask. |
| `--genome PATH` | Reference genome FASTA. |
| `--gtf PATH` | Gene annotation GTF. |
| `--barcode-map PATH` | ST_BarcodeMap executable. |
| `--threads N` | Worker threads; default is 32 or `SLURM_CPUS_PER_TASK`. |
| `--segthreshold X` | Segmentation error threshold; default is 0.15. |
| `--adapter-fa PATH` | Override the packaged adapter FASTA. |
| `--anchor-fa PATH` | Override the packaged anchor FASTA. |
| `--segment-only` | Stop after cDNA segmentation. |
| `--fragment-out PATH` | Explicit segmentation output path. |
| `--skip-index` | Reuse an existing precise CID index. |
| `--skip-fastq` | Reuse an existing converted/decompressed FASTQ. |
| `--run-bambu` | Run Bambu transcript annotation. |
| `--merge-annot` | Merge Bambu assignments with spatial CID results. |
| `--config-env PATH` | Load site-specific defaults. |
| `--version` | Print the toolkit version. |

Run `fullscope run --help` for the complete parser-supported interface,
including explicit intermediate-file overrides.

## Unified command

Fullscope has one public command with two processing modes:

```text
fullscope run [options]       complete spatial-transcriptomics workflow
fullscope segment [options]   cDNA segmentation only
```

For convenience, complete-workflow options can also be passed directly, for
example `fullscope --sample sample01 ...`. The compiled C++ engine is installed
as an internal component and is not part of the public command-line interface.

Starting with version 1.2.0, `fullscope run` replaces `fullscope-ont`, and
`fullscope segment` replaces `fullscope-segment`.

The toolkit uses the precise `build_idx p` and `map_p` path by default.

## Slurm

The packaged Slurm script intentionally contains no cluster-specific partition
or memory setting. Supply those values through `FULLSCOPE_SBATCH_ARGS`:

```bash
FULLSCOPE_SBATCH_ARGS="--partition=compute --mem=300G --cpus-per-task=32" \
fullscope-submit \
  --sample sample01 \
  --outdir results/sample01 \
  --raw-fq reads.fastq.gz \
  --stereoindex sample01.barcodeToPos.h5 \
  --config-env site.env
```

For multiple samples:

```bash
cp fullscope_toolkit/config/samples.example.tsv samples.tsv

fullscope-batch \
  --config-env site.env \
  --samples samples.tsv \
  --sbatch-extra "--partition=compute --mem=300G"
```

## Outputs

A complete run creates:

```text
results/sample01/
  Index/          CID whitelist and precise index
  Fqsegment/      segmented full-length cDNA reads
  CIDextract/     per-read CID candidates
  CIDmap/         mapped spatial CIDs
  Alignment/      splice-aware BAM and index
  Bambu/          optional transcript assignments
  raw_fastq/      BAM-converted or decompressed FASTQ, when needed
  *_fsraw_merged_data.qs
  *_fsraw_merged_data_uniquereads.qs
```

The two `.qs` matrices contain transcript/gene annotation and spatial `x`, `y`
coordinates. They are produced only when `--run-bambu --merge-annot` is used.

## Repository layout

```text
Fullscope/
  scripts/src/                 C++ implementation
  scripts/include/             C++ headers
  fullscope_toolkit/
    bin/                       user-facing command wrappers
    scripts/                   workflow, Slurm and R scripts
    config/                    portable configuration examples
    refdata/                   default adapter and anchor sequences
    testdata/                  small smoke-test FASTQ
  tests/smoke_test.sh
  analysis_script/             study analysis notebooks
  Segmentation_script/         legacy segmentation workflow
  CMakeLists.txt
  environment.yml
  install.sh
```

## Reuse notes

- Fullscope does not bundle reference genomes, gene annotations, Stereo-seq
  mask files or ST_BarcodeMap.
- Memory and runtime depend on read count, CID-whitelist size, thread count and
  the selected long-read platform. Validate a small subset before a full run.
- The bundled smoke data validates installation and segmentation, not a
  biological end-to-end result.
- Preserve the editable configuration used for each run and record the
  Fullscope version with `fullscope --version`.

## Archived version

The code associated with the publication is archived at
[Zenodo](https://doi.org/10.5281/zenodo.19647550).

## License

Fullscope is distributed under the [Apache License 2.0](LICENSE).
