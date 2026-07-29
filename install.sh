#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash install.sh [options]

Options:
  --prefix PATH       Installation prefix.
                      Default: $CONDA_PREFIX when active, otherwise $HOME/.local
  --build-dir PATH    CMake build directory. Default: ./build
  --jobs N            Parallel build jobs. Default: detected CPU count
  -h, --help          Show this help message.
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
prefix="${CONDA_PREFIX:-${HOME}/.local}"
build_dir="${repo_root}/build"
jobs="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix) prefix="$2"; shift 2 ;;
        --build-dir) build_dir="$2"; shift 2 ;;
        --jobs) jobs="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

command -v cmake >/dev/null 2>&1 || {
    echo "ERROR: cmake is required. Create the conda environment first:" >&2
    echo "  conda env create -f environment.yml" >&2
    exit 1
}

generator_args=()
if [[ ! -f "${build_dir}/CMakeCache.txt" ]] && command -v ninja >/dev/null 2>&1; then
    generator_args=(-G Ninja)
fi

cmake -S "$repo_root" -B "$build_dir" "${generator_args[@]}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$prefix"
cmake --build "$build_dir" --parallel "$jobs"
cmake --install "$build_dir"

cat <<EOF
Fullscope installed successfully.

Prefix:
  $prefix

Commands:
  $prefix/bin/fullscope
  $prefix/bin/fullscope-ont
  $prefix/bin/fullscope-segment
  $prefix/bin/fullscope-submit

If needed, add the installation to PATH:
  export PATH="$prefix/bin:\$PATH"
EOF
