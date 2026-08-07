#!/usr/bin/env bash
# =============================================================================
#  phasefield_gpu.sh  --  build and/or run the 3D phase-field solver on an
#  NVIDIA GPU using the NVIDIA HPC SDK (nvfortran) under WSL2 Ubuntu.
#
#  Usage:
#     bash phasefield_gpu.sh              # build (if needed) then run
#     bash phasefield_gpu.sh build        # compile only
#     bash phasefield_gpu.sh run          # run only (must be built already)
#     bash phasefield_gpu.sh rebuild      # force a clean recompile, then run
#     bash phasefield_gpu.sh run 60       # run, but stop after 60 s (quick test)
#
#  Everything is auto-detected: the HPC SDK version, the bundled CUDA toolkit,
#  and the GPU's compute capability. Override any of them with env vars, e.g.
#     CC=80 bash phasefield_gpu.sh        # force compute capability 8.0
#     STRICT=0 bash phasefield_gpu.sh     # allow silent CPU fallback (see below)
# =============================================================================
set -euo pipefail

EXE=phasefield3d_gpu_wsl
SRC="source/variaveis.f90 source/auxiliar.f90 source/concentracao_soluto.f90 \
source/diagnosticos.f90 source/entrada.f90 source/equacao_campo_fase.f90 \
source/principal.f90 source/saida.f90 source/temperatura.f90 source/W0_E0.f90"

# --- run from the directory the script lives in -----------------------------
cd "$(dirname "$0")"

# --- locate the NVIDIA HPC SDK (newest numeric version wins) ----------------
HPC_BASE=/opt/nvidia/hpc_sdk/Linux_x86_64
if [ ! -d "$HPC_BASE" ]; then
    echo "ERROR: NVIDIA HPC SDK not found under $HPC_BASE" >&2
    echo "       Install it, or edit HPC_BASE in this script." >&2
    exit 1
fi
HPC_ROOT=$(ls -d "$HPC_BASE"/*/ 2>/dev/null | grep -E '/[0-9]+\.[0-9]+/$' | sort -V | tail -1)
HPC_ROOT=${HPC_ROOT%/}
NVBIN="$HPC_ROOT/compilers/bin"
if [ ! -x "$NVBIN/nvfortran" ]; then
    echo "ERROR: nvfortran not found in $NVBIN" >&2
    exit 1
fi

# --- pick the CUDA toolkit bundled with that SDK (newest) -------------------
CUDA_VER=$(ls -d "$HPC_ROOT"/cuda/*/ 2>/dev/null | grep -E '/[0-9]+\.[0-9]+/$' | sort -V | tail -1)
CUDA_VER=$(basename "${CUDA_VER%/}")
CUDA_LIB="$HPC_ROOT/cuda/$CUDA_VER/lib64"

# --- detect GPU compute capability (nvidia-smi -> e.g. 8.9 -> cc89) ---------
#     Override with:  CC=89 bash phasefield_gpu.sh
if [ -z "${CC:-}" ]; then
    CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d ' .')
    CC=${CAP:-89}          # fall back to 89 (RTX 40xx) if the query fails
fi

export PATH="$NVBIN:$PATH"
# CRITICAL for WSL2: the CUDA *driver* (libcuda.so.1) lives in /usr/lib/wsl/lib
# and is NOT on the SDK's default search path. Without it the OpenMP runtime
# reports 0 devices and silently falls back to the CPU path.
export LD_LIBRARY_PATH="/usr/lib/wsl/lib:$CUDA_LIB:$HPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"

echo "HPC SDK : $HPC_ROOT"
echo "CUDA    : $CUDA_VER"
echo "GPU     : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1) (cc$CC)"
echo

build() {
    echo ">>> building $EXE (nvfortran -mp=gpu -gpu=cc$CC,cuda$CUDA_VER -O3)"
    rm -f ./*.o ./*.mod
    nvfortran -mp=gpu -gpu="cc$CC,cuda$CUDA_VER" -O3 -Minfo=mp $SRC -o "$EXE"
    rm -f ./*.o ./*.mod
    echo ">>> build ok -> ./$EXE"
    echo
}

run() {
    [ -x "./$EXE" ] || { echo "ERROR: ./$EXE not built yet -- run: bash $0 build" >&2; exit 1; }
    mkdir -p output                      # saida.f90 writes into ./output/
    # STRICT=1 (default): abort loudly if the kernels can't run on the GPU,
    # instead of silently using the CPU. Set STRICT=0 to allow CPU fallback.
    if [ "${STRICT:-1}" = "1" ]; then
        export OMP_TARGET_OFFLOAD=MANDATORY
    fi
    echo ">>> running ./$EXE   (Ctrl-C to stop; output in ./output/)"
    if [ -n "${1:-}" ]; then
        echo ">>> (auto-stop after ${1}s)"
        timeout "$1" "./$EXE"
    else
        "./$EXE"
    fi
}

case "${1:-all}" in
    build)          build ;;
    rebuild)        build; run "${2:-}" ;;
    run)            run   "${2:-}" ;;
    all|"")         [ -x "./$EXE" ] || build; run "${2:-}" ;;
    *) echo "usage: bash $0 [build|run|rebuild|all] [timeout_seconds]" >&2; exit 2 ;;
esac
