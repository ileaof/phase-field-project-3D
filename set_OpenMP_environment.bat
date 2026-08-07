@echo off

rem call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

set OMP_NUM_THREADS=8
set OMP_PROC_BIND=close
set OMP_PLACES=cores

rem Alternative (Intel runtime)
rem set KMP_AFFINITY=granularity=fine,compact,1,0

show_OpenMP_environment.bat