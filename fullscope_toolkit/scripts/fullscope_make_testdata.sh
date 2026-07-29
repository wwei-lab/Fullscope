#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  fullscope_make_testdata.sh --input-fq input.fastq --output-fq test.fastq [options]

Required:
  --input-fq PATH       Source FASTQ
  --output-fq PATH      Downsampled FASTQ

Optional:
  --reads N             Number of reads to keep. Default: 1000
  --version             Print toolkit version and exit

Example:
  fullscope-make-testdata \
    --input-fq /data/longreads/sample.fastq \
    --output-fq /data/testdata/sample_1k.fastq \
    --reads 1000
EOF
}

input_fq=""
output_fq=""
reads=1000

toolkit_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
toolkit_version_file="${toolkit_root}/VERSION"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-fq) input_fq="$2"; shift 2 ;;
        --output-fq) output_fq="$2"; shift 2 ;;
        --reads) reads="$2"; shift 2 ;;
        --version)
            if [[ -f "${toolkit_version_file}" ]]; then
                cat "${toolkit_version_file}"
            else
                echo "unknown"
            fi
            exit 0
            ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ -z "${input_fq}" || -z "${output_fq}" ]]; then
    echo "ERROR: --input-fq and --output-fq are required." >&2
    usage >&2
    exit 1
fi

[[ -e "${input_fq}" ]] || { echo "ERROR: file not found: ${input_fq}" >&2; exit 1; }
[[ "${reads}" =~ ^[0-9]+$ ]] || { echo "ERROR: --reads must be a positive integer." >&2; exit 1; }
if [[ "${reads}" -le 0 ]]; then
    echo "ERROR: --reads must be > 0." >&2
    exit 1
fi

mkdir -p "$(dirname "${output_fq}")"
line_count=$((reads * 4))

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building test FASTQ"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] input=${input_fq}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] output=${output_fq}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] reads=${reads}"

if [[ "$input_fq" == *.gz ]]; then
    command -v gzip >/dev/null 2>&1 || { echo "ERROR: gzip is not available in PATH." >&2; exit 1; }
    gzip -dc -- "$input_fq" | sed -n "1,${line_count}p" > "$output_fq"
else
    sed -n "1,${line_count}p" "$input_fq" > "$output_fq"
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done"
