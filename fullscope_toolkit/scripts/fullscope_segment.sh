#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
toolkit_root="$(cd "${script_dir}/.." && pwd)"
toolkit_version_file="${toolkit_root}/VERSION"

usage() {
    cat <<'EOF'
Usage:
  fullscope_segment.sh --raw-fq input.fastq --out fragment.tsv [options]

Required:
  --raw-fq PATH         Input ONT FASTQ
  --out PATH            Output fragment TSV from process_fq

Optional:
  --threads N           Default: 32
  --segthreshold FLOAT  Default: 0.15
  --adapter-fa PATH     Override built-in adapter FASTA
  --anchor-fa PATH      Override built-in anchor FASTA
  --config-env PATH     Source site/local env config before resolving defaults
  --version             Print toolkit version and exit

Example:
  fullscope segment \
    --raw-fq /data/longreads/sample.fastq \
    --out /data/results/sample_fragment.fastq
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

raw_fq=""
out_tsv=""
threads="${SLURM_CPUS_PER_TASK:-32}"
segthreshold=0.15

adapter_fa="${ADAPTER_FA:-${toolkit_root}/refdata/adapters.fa}"
anchor_fa="${ANCHOR_FA:-${toolkit_root}/refdata/anchor.fa}"
core_bin="${FULLSCOPE_CORE_BIN:-}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --raw-fq) raw_fq="$2"; shift 2 ;;
        --out) out_tsv="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        --segthreshold) segthreshold="$2"; shift 2 ;;
        --adapter-fa) adapter_fa="$2"; shift 2 ;;
        --anchor-fa) anchor_fa="$2"; shift 2 ;;
        --core-bin) core_bin="$2"; shift 2 ;;
        --config-env) config_env="$2"; shift 2 ;;
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

if [[ -z "${raw_fq}" || -z "${out_tsv}" ]]; then
    echo "ERROR: --raw-fq and --out are required." >&2
    usage >&2
    exit 1
fi

[[ -e "${raw_fq}" ]] || { echo "ERROR: file not found: ${raw_fq}" >&2; exit 1; }
[[ -e "${adapter_fa}" ]] || { echo "ERROR: file not found: ${adapter_fa}" >&2; exit 1; }
[[ -e "${anchor_fa}" ]] || { echo "ERROR: file not found: ${anchor_fa}" >&2; exit 1; }
[[ -n "${core_bin}" ]] || { echo "ERROR: the internal Fullscope core is unavailable." >&2; exit 1; }
[[ -x "${core_bin}" ]] || { echo "ERROR: internal Fullscope core is not executable: ${core_bin}" >&2; exit 1; }

mkdir -p "$(dirname "${out_tsv}")"
ulimit -n 65535 || true

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running process_fq"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] raw_fq=${raw_fq}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] out=${out_tsv}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] threads=${threads}"

"${core_bin}" process_fq "${raw_fq}" "${adapter_fa}" "${anchor_fa}" "${segthreshold}" "${threads}" "${out_tsv}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done"
