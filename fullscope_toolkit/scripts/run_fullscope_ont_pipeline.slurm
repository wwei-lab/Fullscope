#!/usr/bin/env bash
#SBATCH -J fullscope_pipe
#SBATCH -N 1
#SBATCH --output=fullscope-%x-%j.out
#SBATCH --error=fullscope-%x-%j.err
#SBATCH --ntasks=1
#SBATCH -c 32

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
installed_core="${SCRIPT_DIR}/../../../libexec/fullscope/fullscope-core"
source_core="${SCRIPT_DIR}/../../build/fullscope-core"
if [[ -x "${installed_core}" ]]; then
    export FULLSCOPE_CORE_BIN="${FULLSCOPE_CORE_BIN:-${installed_core}}"
elif [[ -x "${source_core}" ]]; then
    export FULLSCOPE_CORE_BIN="${FULLSCOPE_CORE_BIN:-${source_core}}"
fi
exec "${SCRIPT_DIR}/fullscope_ont_pipeline.sh" "$@"
