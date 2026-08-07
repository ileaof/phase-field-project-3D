Compiler flags (Windows ifx, matching your existing compile.bat)

ifx -Qopenmp -O2 -QxHost source\variaveis.f90 source\auxiliar.f90 source\concentracao_soluto.f90 ^
    source\diagnosticos.f90 source\entrada.f90 source\equacao_campo_fase.f90 source\principal.f90 ^
    source\saida.f90 source\temperatura.f90 source\W0_E0.f90 -o phasefield3d_omp.exe
(-qopenmp lowercase is rejected by ifx on Windows — needs capital -Qopenmp.) Runtime requires libiomp5md.dll (in the Intel compiler bin folder) on PATH. I did not add -Qipo/-fp:fast=2 — cross-TU inlining and fast-math reordering interact with the OpenMP reduction/RNG determinism claims above, so treat those as an opt-in experiment to re-verify (re-run the thread-count diff test) rather than a default.

Recommended next step for you
Set OMP_PROC_BIND=close and OMP_PLACES=cores (or KMP_AFFINITY=granularity=fine,compact,1,0) before production runs — the 8-thread dip I measured (worse than 4-thread) is a classic symptom of unpinned threads on this machine, and this should resolve it. Verification tools: use diagnosticos.f90's L2/energy/min-max/solute-total output (now printed automatically at every write checkpoint) to compare serial vs. parallel runs beyond what I could test here.

c:\Users\ileao\OneDrive\Documentos\phasefield\build.ps1
create a para parallel compilation file at root dir called parallel_compile.bat
Read c:\Users\ileao\OneDrive\Documentos\phasefield\compile.bat
Read c:\Users\ileao\OneDrive\Documentos\phasefield\build.ps1
Read c:\Users\ileao\OneDrive\Documentos\phasefield\setvars.bat
Write c:\Users\ileao\OneDrive\Documentos\phasefield\parallel_compile.bat
2 lines
Let me verify the batch file actually works standalone:

Bash
IN
cd "c:/Users/ileao/OneDrive/Documentos/phasefield"
export LIB="C:\Program Files (x86)\Intel\oneAPI\compiler\2024.2\lib;C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\ATLMFC\lib\x64;C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\lib\x64;C:\Program Files (x86)\Windows Kits\NETFXSDK\4.8\lib\um\x64;C:\Program Files (x86)\Windows Kits\10\lib\10.0.26100.0\ucrt\x64;C:\Program Files (x86)\Windows Kits\10\lib\10.0.26100.0\um\x64"
export PATH="/c/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.44.35207/bin/Hostx64/x64:/c/Program Files (x86)/Intel/oneAPI/compiler/2024.2/bin:/c/Program Files (x86)/Intel/oneAPI/compiler/2024.2/bin/intel64:$PATH"
rm -f phasefield3d_omp.exe
cmd //c parallel_compile.bat 2>&1 | grep -v "remark #8291\|^write(\|^---"
ls -la phasefield3d_omp.exe

OUT
Exit code 2
'parallel_compile.bat' n o   reconhecido como um comando interno
ou externo, um programa oper vel ou um arquivo em lotes.
ls: cannot access 'phasefield3d_omp.exe': No such file or directory
Bash
IN
cd "c:/Users/ileao/OneDrive/Documentos/phasefield"
export LIB="C:\Program Files (x86)\Intel\oneAPI\compiler\2024.2\lib;C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\ATLMFC\lib\x64;C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\lib\x64;C:\Program Files (x86)\Windows Kits\NETFXSDK\4.8\lib\um\x64;C:\Program Files (x86)\Windows Kits\10\lib\10.0.26100.0\ucrt\x64;C:\Program Files (x86)\Windows Kits\10\lib\10.0.26100.0\um\x64"
export PATH="/c/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.44.35207/bin/Hostx64/x64:/c/Program Files (x86)/Intel/oneAPI/compiler/2024.2/bin:$PATH"
cmd //c 'C:\Users\ileao\OneDrive\Documentos\phasefield\parallel_compile.bat' 2>&1 | grep -v "remark #8291\|^write(\|^---"
ls -la phasefield3d_omp.exe

OUT
C:\Users\ileao\OneDrive\Documentos\phasefield>ifx -Qopenmp -O2 source\variaveis.f90 source\auxiliar.f90 source\concentracao_soluto.f90 source\diagnosticos.f90 source\entrada.f90 source\equacao_campo_fase.f90 source\principal.f90 source\saida.f90 source\temperatura.f90 source\W0_E0.f90 -o phasefield3d_omp.exe 
Intel(R) Fortran Compiler for applications running on Intel(R) 64, Version 2024.2.1 Build 20240711
Copyright (C) 1985-2024 Intel Corporation. All rights reserved.

Microsoft (R) Incremental Linker Version 14.44.35228.0
Copyright (C) Microsoft Corporation.  All rights reserved.

-out:phasefield3d_omp.exe 
-subsystem:console 
-defaultlib:libiomp5md.lib 
-nodefaultlib:vcomp.lib 
-nodefaultlib:vcompd.lib 
variaveis.obj 
auxiliar.obj 
concentracao_soluto.obj 
diagnosticos.obj 
entrada.obj 
equacao_campo_fase.obj 
principal.obj 
saida.obj 
temperatura.obj 
W0_E0.obj 
-rwxr-xr-x 1 ileao 197609 804864 Jul  5 14:04 phasefield3d_omp.exe
parallel_compile.bat compiles cleanly and produces phasefield3d_omp.exe. Cleaning up build artifacts:

Bash
IN
cd "c:/Users/ileao/OneDrive/Documentos/phasefield"
rm -f *.obj *.mod
ls

OUT
build.ps1
compile.bat
entrada.dat
output
parallel_compile.bat
phasefield3d.exe
phasefield3d_omp.exe
setvars.bat
source
Created parallel_compile.bat at the project root, matching compile.bat's style:


ifx -Qopenmp -O2 source\variaveis.f90 source\auxiliar.f90 source\concentracao_soluto.f90 source\diagnosticos.f90 source\entrada.f90 source\equacao_campo_fase.f90 source\principal.f90 source\saida.f90 source\temperatura.f90 source\W0_E0.f90 -o phasefield3d_omp.exe
Verified it runs standalone and produces a working phasefield3d_omp.exe. Same usage as build.ps1/compile.bat — run it after setvars.bat has set up the Intel oneAPI environment (needs -Qopenmp, capital Q, on Windows; lowercase -qopenmp is silently ignored by ifx there). Note the resulting executable needs libiomp5md.dll (in the Intel compiler's bin folder) on PATH at runtime.

