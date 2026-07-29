#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
binary="${1:-${repo_root}/build/fullscope}"

if [[ ! -x "$binary" ]]; then
    echo "ERROR: Fullscope binary is not executable: $binary" >&2
    exit 1
fi

tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/fullscope-smoke.XXXXXX")"
trap 'rm -rf -- "$tmp_dir"' EXIT

input_fq="${repo_root}/fullscope_toolkit/testdata/smoke_20.fastq"
adapter_fa="${repo_root}/fullscope_toolkit/refdata/adapters.fa"
anchor_fa="${repo_root}/fullscope_toolkit/refdata/anchor.fa"
output_tsv="${tmp_dir}/smoke_fragment.tsv"

"$binary" process_fq \
    "$input_fq" \
    "$adapter_fa" \
    "$anchor_fa" \
    0.15 \
    2 \
    "$output_tsv"

[[ -s "$output_tsv" ]] || {
    echo "ERROR: segmentation output is missing or empty" >&2
    exit 1
}

[[ -s "${output_tsv}.summary.tsv" ]] || {
    echo "ERROR: segmentation summary is missing or empty" >&2
    exit 1
}

echo "Smoke test passed: $(wc -l < "$output_tsv") output rows"
