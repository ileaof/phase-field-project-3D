#!/usr/bin/env bash
# =============================================================================
#  karma_gpu.sh -- compila e roda o karma3d na GPU NVIDIA (nvfortran, WSL2).
#  Uso:
#     bash karma_gpu.sh build     # so compila
#     bash karma_gpu.sh run       # so roda (le karma.nml desta pasta)
#     bash karma_gpu.sh           # compila se preciso e roda
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
EXE=karma3d_gpu

HPC_ROOT=$(ls -d /opt/nvidia/hpc_sdk/Linux_x86_64/*/ | grep -E '/[0-9]+\.[0-9]+/$' | sort -V | tail -1)
HPC_ROOT=${HPC_ROOT%/}
CUDA_VER=$(basename "$(ls -d "$HPC_ROOT"/cuda/*/ | grep -E '/[0-9]+\.[0-9]+/$' | sort -V | tail -1)")
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d ' .')
CC=${CC:-${CAP:-89}}
export PATH="$HPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="/usr/lib/wsl/lib:$HPC_ROOT/cuda/$CUDA_VER/lib64:$HPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"

build() {
    nvfortran -mp=gpu -gpu="cc$CC,cuda$CUDA_VER" -O3 -Minfo=mp karma3d.f90 -o "$EXE" 2>&1 | grep -iE "error|warn" || true
    rm -f ./*.mod ./*.o
    [ -x "./$EXE" ] && echo ">>> build ok -> ./$EXE"
}

run() {
    mkdir -p output
    export OMP_TARGET_OFFLOAD=MANDATORY
    "./$EXE"
}

case "${1:-all}" in
    build) build ;;
    run)   run ;;
    all)   [ -x "./$EXE" ] || build; run ;;
    *) echo "uso: bash $0 [build|run]" >&2; exit 2 ;;
esac
