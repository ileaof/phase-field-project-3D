module variaveis

double precision::               Tm, delta_H, E0, W0, M0, D,D_l,D_s, Cp_l,Cp_s,  ANI, DELTA_EPSLON, DELTA_T, T, D_sol, D_liq, D_soluto_C,dif_conc
double precision, allocatable :: phase_field_old(:,:,:),phase_field_new(:,:,:),XC(:,:,:), YC(:,:,:), ZC(:,:,:), TNEW(:,:,:),TOLD(:,:,:), DPP(:,:,:),conc_soluto_new(:,:,:), conc_soluto_old(:,:,:)
double PRECISION, allocatable :: T_(:,:,:)
double precision              :: XX, YY,XL1,YL1,ZL1, h, h_E, h_W, h_N, h_S, h_T, h_B, CSP, CSN, CSS, CSE, CSW, CST, CSB, CLP, CLN, CLS, CLE, CLW, CLT, CLB, Cl_eq, Cs_eq, Vm, R
double precision::               dx, dy, dz, Me, C0,InvDX,InvDY,InvDZ,Hphi
double precision::               g_linha, h_linha, PHI, small, AA, GG
double precision::               PP,PN,PS,PE,PW,PT,PB,PPNEW, AAA,ruidos, FONTE, TNEW1, T0, TE, TP, TW, TS, TN, TT, TB, CPP, CPS, CPN, CPW, CPE, CPT, CPB
double precision::               GAMA_P, GAMA_S, GAMA_N, GAMA_W, GAMA_E, GAMA_T, GAMA_B, GAMAn, GAMAw, GAMAe, GAMAs, GAMAt, GAMAb, c_new, Thin_Interface, Lambda
! ---- cubic (quartic) 3D anisotropy locals: gradient components + eps/cross-terms at the
! six face-neighbor points (E,W,N,S,T,B), evaluated with the same wide/narrow centered
! finite-difference pattern the 2D code used for the angle-based epsilon at each face.
! The diagonal (variable-coefficient Laplacian) and cross (mixed-derivative) parts are
! assembled with the same narrow-difference / corner-weighted pattern as the original 2D
! AA/BB/CC/DD/EE terms (verified by direct numerical cross-check against the 2D formula),
! now with three cross-plane terms (xy,yz,zx) instead of one.
double precision::               K_aniso, Cross_xy, Cross_yz, Cross_zx
double precision::               phix_E, phiy_E, phiz_E, Q_E, S4_E, eps_E, exy_E, eyz_E, ezx_E
double precision::               phix_W, phiy_W, phiz_W, Q_W, S4_W, eps_W, exy_W, eyz_W, ezx_W
double precision::               phix_N, phiy_N, phiz_N, Q_N, S4_N, eps_N, exy_N, eyz_N, ezx_N
double precision::               phix_S, phiy_S, phiz_S, Q_S, S4_S, eps_S, exy_S, eyz_S, ezx_S
double precision::               phix_T, phiy_T, phiz_T, Q_T, S4_T, eps_T, exy_T, eyz_T, ezx_T
double precision::               phix_B, phiy_B, phiz_B, Q_B, S4_B, eps_B, exy_B, eyz_B, ezx_B
double precision::               PNE, PSE, PNW, PSW, PNT, PNB, PST, PSB, PET, PEB, PWT, PWB
integer ::                        i,j,k, Nx, Ny, Nz, it_max,it,E_R,NYF,n
integer ::                       WriteFileD_timeCounter
real            ::               ran, Ke,time, Sigma, Beta, ajuste,teta_1, ran_1, teta_2,ran_2
real            ::               teta_3,ran_3,ran_0,av,D1,D2,TX

double precision             ::               stoptime
character(len=10) ::             text
character(len=14) ::             text2
character(len=38) ::             text3

! Runtime GPU-offload switch. Set once in principal.f90 from omp_get_num_devices():
! .true.  -> compute-heavy kernels run as OpenMP target regions on the accelerator;
! .false. -> the pre-existing CPU OpenMP (!$omp parallel do) paths run unchanged.
! Declared here (module scope) so every kernel subroutine sees the same value without
! having to thread it through argument lists.
logical :: use_gpu = .false.

! declare target (named form, in the module specification part): emit a DEVICE
! version of pseudo_random so it can be CALLed from inside the GPU target regions
! (phase-field and solute noise). It is a self-contained integer/bit-mixing pure
! function -- no I/O, no module state, no host-only intrinsics -- so it is safe and
! reproducible on the accelerator, and race-free when evaluated concurrently.
!$omp declare target(pseudo_random)

contains

! Stateless, deterministic per-cell noise generator replacing "call random(ran)".
! No shared/mutable state -> race-free under OpenMP, and bit-reproducible regardless
! of thread count or scheduling, since each draw is a pure function of (i,j,k,it,salt)
! rather than of call order in a shared stream. "salt" decorrelates independent call
! sites (e.g. phase-field noise vs concentration noise) that share the same (i,j,k,it).
! Hash core is the MurmurHash3 64-bit finalizer (fmix64), combined with prime-multiplier
! mixing of the four integer inputs -- standard, well-tested bit-mixing constants.
pure function pseudo_random(i,j,k,it,salt) result(ran)
implicit none
integer, intent(in) :: i,j,k,it,salt
double precision :: ran
integer(kind=8) :: hval
integer(kind=8), parameter :: C1 = int(z'ff51afd7ed558ccd',8)
integer(kind=8), parameter :: C2 = int(z'c4ceb9fe1a85ec53',8)
integer(kind=8), parameter :: MASK63 = int(z'7fffffffffffffff',8)

hval = int(i,8)*2654435761_8 + int(j,8)*2246822519_8 + int(k,8)*3266489917_8 &
     + int(it,8)*668265263_8 + int(salt,8)*374761393_8 + 1_8
hval = ieor(hval, ishft(hval,-33))
hval = hval * C1
hval = ieor(hval, ishft(hval,-33))
hval = hval * C2
hval = ieor(hval, ishft(hval,-33))

ran = real(iand(hval, MASK63), 8) / real(MASK63, 8)

end function pseudo_random

end module variaveis




  




