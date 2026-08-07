@echo off
REM ===========================================================================
REM  GPU-offload build (Intel ifx, OpenMP target -> Intel GPU / SPIR64).
REM  Run AFTER the Intel oneAPI environment is set (setvars.bat), exactly like
REM  compile.bat / parallel_compile.bat.
REM
REM  Same single source tree as the CPU build -- the OpenMP `target` regions
REM  fall back to the host when no device is present, and principal.f90 chooses
REM  the CPU or GPU code path at RUNTIME via omp_get_num_devices(). So this exe
REM  runs correctly with OR without a GPU; it just additionally offloads when an
REM  OpenMP target device exists.
REM
REM  -Qiopenmp                 : Intel LLVM-based OpenMP (required for offload)
REM  -Qopenmp-targets:spir64   : generate device code for Intel GPUs
REM  Verify at runtime with:   set LIBOMPTARGET_DEBUG=1   (prints device usage)
REM ===========================================================================
ifx -Qiopenmp -Qopenmp-targets:spir64 -O2 ^
    source\variaveis.f90 source\auxiliar.f90 source\concentracao_soluto.f90 ^
    source\diagnosticos.f90 source\entrada.f90 source\equacao_campo_fase.f90 ^
    source\principal.f90 source\saida.f90 source\temperatura.f90 source\W0_E0.f90 ^
    -o phasefield3d_gpu.exe
