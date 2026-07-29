#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  fullscope_batch_submit.sh --samples samples.tsv [--config-env site.env]

Required:
  --samples FILE         Tab-delimited sample sheet.

Optional:
  --config-env FILE      Source this env file and pass it through to each job.
  --sbatch-extra "..."   Extra sbatch flags prepended before the slurm script.

Sample sheet columns:
  sample  outdir  stereoindex  input_bam  raw_fq  out_prefix  run_bambu  merge_annot  threads  genome  gtf  segment_only  fragment_out

Notes:
  - Keep the header line.
  - Use NA for unused optional values.
  - input_bam or raw_fq must be provided.
EOF
}

samples=""
config_env=""
sbatch_extra=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --samples) samples="$2"; shift 2 ;;
        --config-env) config_env="$2"; shift 2 ;;
        --sbatch-extra) sbatch_extra="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ -z "${samples}" ]]; then
    echo "ERROR: --samples is required." >&2
    exit 1
fi
if [[ ! -s "${samples}" ]]; then
    echo "ERROR: sample sheet not found or empty: ${samples}" >&2
    exit 1
fi
if [[ -n "${config_env}" && ! -f "${config_env}" ]]; then
    echo "ERROR: config env not found: ${config_env}" >&2
    exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
slurm_script="${script_dir}/run_fullscope_ont_pipeline.slurm"

tail -n +2 "${samples}" | while IFS=$'\t' read -r sample outdir stereoindex input_bam raw_fq out_prefix run_bambu merge_annot threads genome gtf segment_only fragment_out; do
    [[ -z "${sample:-}" || "${sample:0:1}" == "#" ]] && continue

    cmd=(sbatch)
    if [[ -n "${sbatch_extra}" ]]; then
        # shellcheck disable=SC2206
        extra_parts=(${sbatch_extra})
        cmd+=("${extra_parts[@]}")
    fi

    cmd+=("${slurm_script}"
        --sample "${sample}"
        --outdir "${outdir}"
        --stereoindex "${stereoindex}")

    if [[ -n "${config_env}" ]]; then
        cmd+=(--config-env "${config_env}")
    fi
    if [[ -n "${input_bam}" && "${input_bam}" != "NA" ]]; then
        cmd+=(--input-bam "${input_bam}")
    fi
    if [[ -n "${raw_fq}" && "${raw_fq}" != "NA" ]]; then
        cmd+=(--raw-fq "${raw_fq}")
    fi
    if [[ -n "${out_prefix}" && "${out_prefix}" != "NA" ]]; then
        cmd+=(--out-prefix "${out_prefix}")
    fi
    if [[ -n "${threads}" && "${threads}" != "NA" ]]; then
        cmd+=(--threads "${threads}")
    fi
    if [[ -n "${genome}" && "${genome}" != "NA" ]]; then
        cmd+=(--genome "${genome}")
    fi
    if [[ -n "${gtf}" && "${gtf}" != "NA" ]]; then
        cmd+=(--gtf "${gtf}")
    fi
    if [[ -n "${fragment_out}" && "${fragment_out}" != "NA" ]]; then
        cmd+=(--fragment-out "${fragment_out}")
    fi
    if [[ "${run_bambu:-0}" == "1" ]]; then
        cmd+=(--run-bambu)
    fi
    if [[ "${merge_annot:-0}" == "1" ]]; then
        cmd+=(--merge-annot)
    fi
    if [[ "${segment_only:-0}" == "1" ]]; then
        cmd+=(--segment-only)
    fi

    echo "Submitting ${sample}"
    "${cmd[@]}"
done
