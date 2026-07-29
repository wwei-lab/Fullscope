#!/usr/bin/env bash
#SBATCH -J fullscope_pipe
#SBATCH -N 1
#SBATCH --output=fullscope-%x-%j.out
#SBATCH --error=fullscope-%x-%j.err
#SBATCH --ntasks=1
#SBATCH -c 32

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${SCRIPT_DIR}/fullscope_ont_pipeline.sh" "$@"
