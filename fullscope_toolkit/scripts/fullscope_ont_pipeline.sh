#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
toolkit_root="$(cd "${script_dir}/.." && pwd)"
toolkit_version_file="${toolkit_root}/VERSION"

usage() {
    cat <<'EOF'
Usage:
  fullscope_ont_pipeline.sh \
    --sample E13.5_mid400 \
    --outdir /path/to/processed_data/E13.5_mid400 \
    --stereoindex /path/to/Y01460J7.barcodeToPos.h5 \
    --input-bam /path/to/input.bam \
    [--run-bambu] [--merge-annot]

Optional:
  --config-env PATH     Source a site/local env file before resolving defaults.
  --segment-only        Run only process_fq and stop after segmentation output.
  --fragment-out PATH   Explicit output path for process_fq TSV.
  --bambu-min-read-length N
                        Minimum query length retained for Bambu (default: 200).
  --bambu-max-read-length N
                        Maximum query length retained for Bambu (default: 20000).
  --bambu-bam PATH      Explicit filtered BAM checkpoint/output path.
  --skip-bambu-filter   Pass the complete alignment directly to Bambu.
  --version             Print toolkit version and exit.

External requirement for CID index creation:
  ST_BarcodeMap must be in PATH or supplied with --barcode-map.
EOF
}

config_env=""
cli_args=("$@")
for ((i = 0; i < ${#cli_args[@]}; i++)); do
    if [[ "${cli_args[$i]}" == "--config-env" ]]; then
        ((i + 1 < ${#cli_args[@]})) || {
            echo "ERROR: --config-env requires a path." >&2
            exit 2
        }
        config_env="${cli_args[$((i + 1))]}"
        break
    fi
done

if [[ -n "${config_env}" ]]; then
    if [[ ! -f "${config_env}" ]]; then
        echo "ERROR: --config-env file not found: ${config_env}" >&2
        exit 1
    fi
    # shellcheck disable=SC1090
    source "${config_env}"
fi

sample=""
outdir=""
stereoindex=""
input_bam=""
raw_fq=""
out_prefix=""

run_bambu=0
merge_annot=0
skip_index=0
skip_fastq=0
segment_only=0
skip_bambu_filter=0

threads="${SLURM_CPUS_PER_TASK:-32}"
kmer=7
bucketnum=6
segthreshold=0.15
bambu_min_read_length=200
bambu_max_read_length=20000

adapter_fa="${ADAPTER_FA:-${toolkit_root}/refdata/adapters.fa}"
anchor_fa="${ANCHOR_FA:-${toolkit_root}/refdata/anchor.fa}"
genome="${GENOME_FA:-}"
gtf="${GENES_GTF:-}"
barcode_map="${ST_BARCODES_MAP:-$(command -v ST_BarcodeMap-0.0.1 || command -v ST_BarcodeMap || true)}"
core_bin="${FULLSCOPE_CORE_BIN:-}"
bambu_r="${BAMBU_R:-${toolkit_root}/scripts/bambu_process.R}"
merge_r="${MERGE_R:-${toolkit_root}/scripts/fs_merge_spatial_transcripts.R}"
rscript_bin="${RSCRIPT_BIN:-$(command -v Rscript || true)}"

cidindex_txt=""
cidindex_prefix=""
cidmap=""
fqalign=""
bambu_bam=""
bambu_out_prefix=""
bambu_sample=""
transdf_qs=""
merged_all_qs=""
merged_unique_qs=""
fragment_out=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config-env) config_env="$2"; shift 2 ;;
        --version)
            if [[ -f "${toolkit_version_file}" ]]; then
                cat "${toolkit_version_file}"
            else
                echo "unknown"
            fi
            exit 0
            ;;
        --segment-only) segment_only=1; shift ;;
        --fragment-out) fragment_out="$2"; shift 2 ;;
        --sample) sample="$2"; shift 2 ;;
        --outdir) outdir="$2"; shift 2 ;;
        --stereoindex) stereoindex="$2"; shift 2 ;;
        --input-bam) input_bam="$2"; shift 2 ;;
        --raw-fq) raw_fq="$2"; shift 2 ;;
        --out-prefix) out_prefix="$2"; shift 2 ;;
        --run-bambu) run_bambu=1; shift ;;
        --merge-annot) merge_annot=1; shift ;;
        --skip-index) skip_index=1; shift ;;
        --skip-fastq) skip_fastq=1; shift ;;
        --threads) threads="$2"; shift 2 ;;
        --kmer) kmer="$2"; shift 2 ;;
        --bucketnum) bucketnum="$2"; shift 2 ;;
        --segthreshold) segthreshold="$2"; shift 2 ;;
        --adapter-fa) adapter_fa="$2"; shift 2 ;;
        --anchor-fa) anchor_fa="$2"; shift 2 ;;
        --genome) genome="$2"; shift 2 ;;
        --gtf) gtf="$2"; shift 2 ;;
        --barcode-map) barcode_map="$2"; shift 2 ;;
        --core-bin) core_bin="$2"; shift 2 ;;
        --bambu-r) bambu_r="$2"; shift 2 ;;
        --bambu-min-read-length) bambu_min_read_length="$2"; shift 2 ;;
        --bambu-max-read-length) bambu_max_read_length="$2"; shift 2 ;;
        --bambu-bam) bambu_bam="$2"; shift 2 ;;
        --skip-bambu-filter) skip_bambu_filter=1; shift ;;
        --merge-r) merge_r="$2"; shift 2 ;;
        --rscript-bin) rscript_bin="$2"; shift 2 ;;
        --cidindex-txt) cidindex_txt="$2"; shift 2 ;;
        --cidindex-prefix) cidindex_prefix="$2"; shift 2 ;;
        --cidmap) cidmap="$2"; shift 2 ;;
        --fqalign) fqalign="$2"; shift 2 ;;
        --bambu-out-prefix) bambu_out_prefix="$2"; shift 2 ;;
        --bambu-sample) bambu_sample="$2"; shift 2 ;;
        --transdf-qs) transdf_qs="$2"; shift 2 ;;
        --merged-all-qs) merged_all_qs="$2"; shift 2 ;;
        --merged-unique-qs) merged_unique_qs="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ -z "$sample" || -z "$outdir" ]]; then
    echo "ERROR: --sample and --outdir are required." >&2
    exit 1
fi
if [[ "$segment_only" -eq 0 && -z "$stereoindex" ]]; then
    echo "ERROR: --sample, --outdir, and --stereoindex are required." >&2
    exit 1
fi
if [[ -z "$input_bam" && -z "$raw_fq" ]]; then
    echo "ERROR: provide either --input-bam or --raw-fq." >&2
    exit 1
fi
if [[ ! "$bambu_min_read_length" =~ ^[0-9]+$ ||
      ! "$bambu_max_read_length" =~ ^[0-9]+$ ||
      "$bambu_min_read_length" -gt "$bambu_max_read_length" ]]; then
    echo "ERROR: Bambu read-length limits must be non-negative integers with min <= max." >&2
    exit 1
fi

out_prefix="${out_prefix:-$sample}"
bambu_sample="${bambu_sample:-$sample}"

mkdir -p "${outdir}/Index" "${outdir}/CIDmap" "${outdir}/CIDextract" "${outdir}/Alignment" "${outdir}/Fqsegment" "${outdir}/raw_fastq" "${outdir}/Bambu"

raw_fq_from_bam=0
if [[ -z "$raw_fq" ]]; then
    raw_fq="${outdir}/raw_fastq/${out_prefix}.fastq"
    raw_fq_from_bam=1
fi
cidindex_txt="${cidindex_txt:-${outdir}/Index/$(basename "${stereoindex}").txt}"
stereo_base="$(basename "${stereoindex}")"
stereo_base="${stereo_base%.h5}"
stereo_base="${stereo_base%.txt}"
cidindex_prefix="${cidindex_prefix:-${outdir}/Index/${stereo_base}_revo}"
freg_fq="${fragment_out:-${outdir}/Fqsegment/${out_prefix}_fragment.tsv}"
cidextract="${outdir}/CIDextract/${out_prefix}_ont_cidextract.tsv"
cidmap="${cidmap:-${outdir}/CIDmap/${out_prefix}}"
fqalign="${fqalign:-${outdir}/Alignment/${out_prefix}_ont.sorted.bam}"
bambu_bam="${bambu_bam:-${fqalign%.bam}.bambu.filtered.q${bambu_min_read_length}-${bambu_max_read_length}.bam}"
bambu_out_prefix="${bambu_out_prefix:-${outdir}/Bambu/${sample}}"
transdf_qs="${transdf_qs:-${bambu_out_prefix}_trans_total_anno.qs}"
merged_all_qs="${merged_all_qs:-${outdir}/${sample}_fsraw_merged_data.qs}"
merged_unique_qs="${merged_unique_qs:-${outdir}/${sample}_fsraw_merged_data_uniquereads.qs}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command $1" >&2; exit 1; }; }
need_file() { [[ -e "$1" ]] || { echo "ERROR: file not found: $1" >&2; exit 1; }; }

ulimit -n 65535 || true
if [[ "$raw_fq_from_bam" -eq 1 ]]; then
    need_cmd bedtools
fi
if [[ "$segment_only" -eq 0 ]]; then
    need_cmd minimap2
    need_cmd samtools
    need_file "$stereoindex"
    [[ -n "$genome" ]] || { echo "ERROR: --genome or GENOME_FA is required." >&2; exit 1; }
    [[ -n "$gtf" ]] || { echo "ERROR: --gtf or GENES_GTF is required." >&2; exit 1; }
fi
need_file "$adapter_fa"
need_file "$anchor_fa"
[[ -n "$core_bin" ]] || { echo "ERROR: the internal Fullscope core is unavailable." >&2; exit 1; }
if [[ "$segment_only" -eq 0 ]]; then
    need_file "$genome"
    need_file "$gtf"
    if [[ "$skip_index" -eq 0 || ! -s "${cidindex_prefix}.precise.bin" ]]; then
        [[ -n "$barcode_map" ]] || {
            echo "ERROR: ST_BarcodeMap is required to build a CID index; use --barcode-map." >&2
            exit 1
        }
        need_file "$barcode_map"
    fi
fi
[[ -x "$core_bin" ]] || { echo "ERROR: internal Fullscope core is not executable: $core_bin" >&2; exit 1; }
if [[ "$run_bambu" -eq 1 ]]; then
    need_file "$bambu_r"
    [[ -n "$rscript_bin" ]] || { echo "ERROR: Rscript is required for --run-bambu." >&2; exit 1; }
    need_file "$rscript_bin"
fi
if [[ "$merge_annot" -eq 1 ]]; then
    need_file "$merge_r"
    [[ -n "$rscript_bin" ]] || { echo "ERROR: Rscript is required for --merge-annot." >&2; exit 1; }
    need_file "$rscript_bin"
fi

if [[ "$segment_only" -eq 0 ]]; then
    if [[ "$skip_index" -eq 0 || ! -s "${cidindex_prefix}.precise.bin" ]]; then
        log "Building CID index"
        "$barcode_map" --in "$stereoindex" --out "$cidindex_txt" --action 3 -w "$threads"
        "$core_bin" build_idx p "$cidindex_txt" "$threads" "$kmer" "$bucketnum" "$cidindex_prefix"
    fi
fi

if [[ "$skip_fastq" -eq 0 ]]; then
    if [[ ! -s "$raw_fq" ]]; then
        need_file "$input_bam"
        log "Converting BAM to FASTQ"
        bedtools bamtofastq -i "$input_bam" -fq "$raw_fq"
    fi
fi

need_file "$raw_fq"

if [[ "$raw_fq" == *.gz ]]; then
    need_cmd gzip
    decompressed_fq="${outdir}/raw_fastq/${out_prefix}.fastq"
    if [[ "$skip_fastq" -eq 0 || ! -s "$decompressed_fq" ]]; then
        log "Decompressing FASTQ"
        gzip -dc -- "$raw_fq" > "$decompressed_fq"
    fi
    raw_fq="$decompressed_fq"
fi

log "Running process_fq"
"$core_bin" process_fq "$raw_fq" "$adapter_fa" "$anchor_fa" "$segthreshold" "$threads" "$freg_fq"

if [[ "$segment_only" -eq 1 ]]; then
    log "Segment-only mode finished"
    log "Fragment TSV: $freg_fq"
    exit 0
fi

log "Running extract_fq"
"$core_bin" extract_fq "$freg_fq" "$cidextract" "$threads"
log "Running map_p"
"$core_bin" map_p "$cidextract" "${cidindex_prefix}.precise.bin" "$cidmap" "$threads" "$kmer"
log "Running minimap2 alignment"
minimap2 -K500m -t "$threads" --secondary=no -a -x splice --splice-flank=yes "$genome" "$freg_fq" | samtools sort -@ "$threads" -o "$fqalign"
samtools index "$fqalign"

if [[ "$run_bambu" -eq 1 ]]; then
    if [[ "$skip_bambu_filter" -eq 1 ]]; then
        bambu_input_bam="$fqalign"
        log "Skipping the Bambu alignment filter by request"
    else
        bambu_input_bam="$bambu_bam"
        if [[ -s "$bambu_input_bam" && -s "${bambu_input_bam}.bai" ]] &&
           samtools quickcheck "$bambu_input_bam"; then
            log "Reusing Bambu-filtered alignment: $bambu_input_bam"
        else
            tmp_bambu_bam="${bambu_input_bam}.tmp.$$.bam"
            filter_counts="${bambu_input_bam}.filter_counts.json"
            tmp_filter_counts="${filter_counts}.tmp.$$"
            rm -f -- "$tmp_bambu_bam" "${tmp_bambu_bam}.bai" "$tmp_filter_counts"
            log "Filtering alignment for Bambu: primary mapped reads, query length ${bambu_min_read_length}-${bambu_max_read_length}"
            samtools view -@ "$threads" -b -F 0x804 \
                -e "qlen >= ${bambu_min_read_length} && qlen <= ${bambu_max_read_length}" \
                --save-counts "$tmp_filter_counts" \
                -o "$tmp_bambu_bam" "$fqalign"
            samtools quickcheck "$tmp_bambu_bam"
            samtools index -@ "$threads" "$tmp_bambu_bam"
            mv -- "$tmp_bambu_bam" "$bambu_input_bam"
            mv -- "${tmp_bambu_bam}.bai" "${bambu_input_bam}.bai"
            mv -- "$tmp_filter_counts" "$filter_counts"
            log "Bambu alignment filter counts: $(tr -d '\n' < "$filter_counts")"
        fi
    fi
    log "Running bambu"
    "$rscript_bin" "$bambu_r" "$gtf" "$genome" "$bambu_input_bam" "$bambu_out_prefix" "$bambu_sample"
fi

if [[ "$merge_annot" -eq 1 ]]; then
    need_file "$transdf_qs"
    need_file "$cidmap"
    log "Merging spatial and transcript annotations"
    "$rscript_bin" "$merge_r" "$transdf_qs" "$cidmap" "$merged_all_qs" "$merged_unique_qs" "$gtf"
fi

log "Done"
