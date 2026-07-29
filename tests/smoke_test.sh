#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
command_path="${1:-${repo_root}/fullscope_toolkit/bin/fullscope}"

if [[ ! -x "$command_path" ]]; then
    echo "ERROR: Fullscope command is not executable: $command_path" >&2
    exit 1
fi

tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/fullscope-smoke.XXXXXX")"
trap 'rm -rf -- "$tmp_dir"' EXIT

input_fq="${repo_root}/fullscope_toolkit/testdata/smoke_20.fastq"
output_tsv="${tmp_dir}/smoke_fragment.tsv"

"$command_path" segment \
    --raw-fq "$input_fq" \
    --out "$output_tsv" \
    --segthreshold 0.15 \
    --threads 2

[[ -s "$output_tsv" ]] || {
    echo "ERROR: segmentation output is missing or empty" >&2
    exit 1
}

[[ -s "${output_tsv}.summary.tsv" ]] || {
    echo "ERROR: segmentation summary is missing or empty" >&2
    exit 1
}

echo "Smoke test passed: $(wc -l < "$output_tsv") output rows"
