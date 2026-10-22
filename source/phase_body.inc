! ============================================================================
! phase_body.inc -- per-cell body of the phase-field update (equacao_campo_fase)
! ----------------------------------------------------------------------------
! Single source of truth for the interior update. It is INCLUDEd verbatim from
! BOTH the CPU worksharing loop (!$omp parallel do) and the GPU offload loop
! (!$omp target teams distribute parallel do) in equacao_campo_fase.f90, so the
! two execution paths are guaranteed bit-for-bit identical in their arithmetic
! (only the surrounding directive differs). Do not add persistent state here:
! every name written below is either a loop index or a module-global scalar that
! is listed PRIVATE on both enclosing directives, so each (i,j,k) iteration owns
! its own copy -- race-free on CPU threads and on GPU teams/threads alike.
! ============================================================================

PP = phase_field_old(i,j,k)

PN= phase_field_old(i,j+1,k)

PS= phase_field_old(i,j-1,k)

PE= phase_field_old(i+1,j,k)

PW= phase_field_old(i-1,j,k)

PT= phase_field_old(i,j,k+1)

PB= phase_field_old(i,j,k-1)

PNE = phase_field_old(i+1,j+1,k)
PSE = phase_field_old(i+1,j-1,k)
PNW = phase_field_old(i-1,j+1,k)
PSW = phase_field_old(i-1,j-1,k)

PNT = phase_field_old(i,j+1,k+1)
PNB = phase_field_old(i,j+1,k-1)
PST = phase_field_old(i,j-1,k+1)
PSB = phase_field_old(i,j-1,k-1)

PET = phase_field_old(i+1,j,k+1)
PEB = phase_field_old(i+1,j,k-1)
PWT = phase_field_old(i-1,j,k+1)
PWB = phase_field_old(i-1,j,k-1)

PPNEW= phase_field_new(i,j,k)

!TP = Told(i,j,k)

small =1.d-23

K_aniso = 16.d0*E0*DELTA_EPSLON

!==============================================PONTO E (face i+1/2,j,k)===============================================
! Local gradient AT node (i+1,j,k): wide centered difference in x (own direction),
! standard centered differences in y,z (mirrors the 2D code's PHIx_E/PHIy_E pattern).

phix_E = (phase_field_old(i+2,j,k)-PP)*0.5d0/dx
phiy_E = (phase_field_old(i+1,j+1,k)-phase_field_old(i+1,j-1,k))*0.5d0/dy
phiz_E = (phase_field_old(i+1,j,k+1)-phase_field_old(i+1,j,k-1))*0.5d0/dz

Q_E  = phix_E**2 + phiy_E**2 + phiz_E**2 + small
S4_E = phix_E**4 + phiy_E**4 + phiz_E**4

eps_E = E0*(1.d0 + DELTA_EPSLON*(4.d0*S4_E/(Q_E**2) - 3.d0))
exy_E = K_aniso*eps_E*phix_E*phiy_E*(phix_E**2-phiy_E**2)/(Q_E**2)
ezx_E = K_aniso*eps_E*phiz_E*phix_E*(phiz_E**2-phix_E**2)/(Q_E**2)

!=============================================PONTO W (face i-1/2,j,k)================================================

phix_W = (PP-phase_field_old(i-2,j,k))*0.5d0/dx
phiy_W = (phase_field_old(i-1,j+1,k)-phase_field_old(i-1,j-1,k))*0.5d0/dy
phiz_W = (phase_field_old(i-1,j,k+1)-phase_field_old(i-1,j,k-1))*0.5d0/dz

Q_W  = phix_W**2 + phiy_W**2 + phiz_W**2 + small
S4_W = phix_W**4 + phiy_W**4 + phiz_W**4

eps_W = E0*(1.d0 + DELTA_EPSLON*(4.d0*S4_W/(Q_W**2) - 3.d0))
exy_W = K_aniso*eps_W*phix_W*phiy_W*(phix_W**2-phiy_W**2)/(Q_W**2)
ezx_W = K_aniso*eps_W*phiz_W*phix_W*(phiz_W**2-phix_W**2)/(Q_W**2)

!============================================PONTO N (face i,j+1/2,k)=================================================

phix_N = (phase_field_old(i+1,j+1,k)-phase_field_old(i-1,j+1,k))*0.5d0/dx
phiy_N = (phase_field_old(i,j+2,k)-PP)*0.5d0/dy
phiz_N = (phase_field_old(i,j+1,k+1)-phase_field_old(i,j+1,k-1))*0.5d0/dz

Q_N  = phix_N**2 + phiy_N**2 + phiz_N**2 + small
S4_N = phix_N**4 + phiy_N**4 + phiz_N**4

eps_N = E0*(1.d0 + DELTA_EPSLON*(4.d0*S4_N/(Q_N**2) - 3.d0))
exy_N = K_aniso*eps_N*phix_N*phiy_N*(phix_N**2-phiy_N**2)/(Q_N**2)
eyz_N = K_aniso*eps_N*phiy_N*phiz_N*(phiy_N**2-phiz_N**2)/(Q_N**2)

!============================================PONTO S (face i,j-1/2,k)=================================================

phix_S = (phase_field_old(i+1,j-1,k)-phase_field_old(i-1,j-1,k))*0.5d0/dx
phiy_S = (PP-phase_field_old(i,j-2,k))*0.5d0/dy
phiz_S = (phase_field_old(i,j-1,k+1)-phase_field_old(i,j-1,k-1))*0.5d0/dz

Q_S  = phix_S**2 + phiy_S**2 + phiz_S**2 + small
S4_S = phix_S**4 + phiy_S**4 + phiz_S**4

eps_S = E0*(1.d0 + DELTA_EPSLON*(4.d0*S4_S/(Q_S**2) - 3.d0))
exy_S = K_aniso*eps_S*phix_S*phiy_S*(phix_S**2-phiy_S**2)/(Q_S**2)
eyz_S = K_aniso*eps_S*phiy_S*phiz_S*(phiy_S**2-phiz_S**2)/(Q_S**2)

!============================================PONTO T (face i,j,k+1/2)=================================================

phix_T = (phase_field_old(i+1,j,k+1)-phase_field_old(i-1,j,k+1))*0.5d0/dx
phiy_T = (phase_field_old(i,j+1,k+1)-phase_field_old(i,j-1,k+1))*0.5d0/dy
phiz_T = (phase_field_old(i,j,k+2)-PP)*0.5d0/dz

Q_T  = phix_T**2 + phiy_T**2 + phiz_T**2 + small
S4_T = phix_T**4 + phiy_T**4 + phiz_T**4

eps_T = E0*(1.d0 + DELTA_EPSLON*(4.d0*S4_T/(Q_T**2) - 3.d0))
eyz_T = K_aniso*eps_T*phiy_T*phiz_T*(phiy_T**2-phiz_T**2)/(Q_T**2)
ezx_T = K_aniso*eps_T*phiz_T*phix_T*(phiz_T**2-phix_T**2)/(Q_T**2)

!============================================PONTO B (face i,j,k-1/2)=================================================

phix_B = (phase_field_old(i+1,j,k-1)-phase_field_old(i-1,j,k-1))*0.5d0/dx
phiy_B = (phase_field_old(i,j+1,k-1)-phase_field_old(i,j-1,k-1))*0.5d0/dy
phiz_B = (PP-phase_field_old(i,j,k-2))*0.5d0/dz

Q_B  = phix_B**2 + phiy_B**2 + phiz_B**2 + small
S4_B = phix_B**4 + phiy_B**4 + phiz_B**4

eps_B = E0*(1.d0 + DELTA_EPSLON*(4.d0*S4_B/(Q_B**2) - 3.d0))
eyz_B = K_aniso*eps_B*phiy_B*phiz_B*(phiy_B**2-phiz_B**2)/(Q_B**2)
ezx_B = K_aniso*eps_B*phiz_B*phix_B*(phiz_B**2-phix_B**2)/(Q_B**2)

!=========================================EQUA��O DE FASE (anisotropic term)=================================

AA = (eps_E**2*(PE-PP)+eps_W**2*(PW-PP))/dx**2 &
   + (eps_N**2*(PN-PP)+eps_S**2*(PS-PP))/dy**2 &
   + (eps_T**2*(PT-PP)+eps_B**2*(PB-PP))/dz**2

Cross_xy = (0.25d0/(dx*dy))*(PNE*(-exy_N-(-exy_E))+PSE*(-exy_E-(-exy_S))+PNW*(-exy_W-(-exy_N))+PSW*(-exy_S-(-exy_W))) &
         + PN*(0.25d0/(dx*dy))*(-exy_W-(-exy_E)) + PS*(0.25d0/(dx*dy))*(-exy_E-(-exy_W)) &
         + PE*(0.25d0/(dx*dy))*(-exy_N-(-exy_S)) + PW*(0.25d0/(dx*dy))*(-exy_S-(-exy_N))

Cross_yz = (0.25d0/(dy*dz))*(PNT*(-eyz_T-(-eyz_N))+PNB*(-eyz_N-(-eyz_B))+PST*(-eyz_S-(-eyz_T))+PSB*(-eyz_B-(-eyz_S))) &
         + PT*(0.25d0/(dy*dz))*(-eyz_S-(-eyz_N)) + PB*(0.25d0/(dy*dz))*(-eyz_N-(-eyz_S)) &
         + PN*(0.25d0/(dy*dz))*(-eyz_T-(-eyz_B)) + PS*(0.25d0/(dy*dz))*(-eyz_B-(-eyz_T))

Cross_zx = (0.25d0/(dz*dx))*(PET*(-ezx_E-(-ezx_T))+PWT*(-ezx_T-(-ezx_W))+PEB*(-ezx_B-(-ezx_E))+PWB*(-ezx_W-(-ezx_B))) &
         + PE*(0.25d0/(dz*dx))*(-ezx_B-(-ezx_T)) + PW*(0.25d0/(dz*dx))*(-ezx_T-(-ezx_B)) &
         + PT*(0.25d0/(dz*dx))*(-ezx_E-(-ezx_W)) + PB*(0.25d0/(dz*dx))*(-ezx_W-(-ezx_E))

AA = AA + Cross_xy + Cross_yz + Cross_zx

!====================================================================================================

g_linha = 2.d0*PP*(1.d0-PP)*(1.d0-PP) - 2.d0*PP*PP*(1.d0-PP)                     !PP*(2-(6*PP)+(4*(PP**2)))

! h(phi) quintico (WBM): h' = 30 phi^2 (1-phi)^2 se anula QUADRATICAMENTE em phi=0,
! entao o liquido (phi=0) e sempre metaestavel. Com o cubico anterior,
! h' = 6 phi (1-phi) e o liquido fica linearmente instavel quando
! 3 (RT/Vm) ln(...) > w0 -- o dominio inteiro solidificava sozinho.
!h_linha = 2.d0*PP*(3.d0-2.d0*PP)-2.d0*PP*PP   ! anterior (cubico)
h_linha = 30.d0*PP*PP*(1.d0-PP)*(1.d0-PP)

CPP  = conc_soluto_old(i,j,k)

!h   = PP*PP*(3.d0-2.d0*PP)                   ! anterior (cubico)
h   = PP*PP*PP*(10.d0-15.d0*PP+6.d0*PP*PP)

CLP  = CPP/(Ke*h+(1.d0-h))

CSP  = Ke*CPP/(Ke*h+(1.d0-h))

! Cs_eq/Cl_eq computed once before the parallel/target region (loop-invariant) -- not
! recomputed per-point here.

GG  = (-w0*g_linha )+ (R*TP/Vm)*h_linha*Log(abs(((1.d0-CS_eq)*(1.d0-CLP))/((1.d0-CL_eq)*(1.d0-CSP))))


!========================================Geracao de ruidos===========================================

 !ran = pseudo_random(i,j,k,it,1)   ! anterior: ruido celula a celula (enruga a interface)
 ran = ruido_suave(i,j,k,it)

 ! Ruido de MEDIA ZERO: pseudo_random devolve [0,1); com ran direto o ruido era
 ! sempre positivo e agia como forca motriz isotropica extra (crescimento esferico).
 !ruidos= ran*AAA*16.d0*(PP**2)*(1-PP)**2.d0   ! anterior (media positiva)
 ruidos= (2.d0*ran-1.d0)*AAA*16.d0*(PP**2)*(1-PP)**2.d0



PPNEW = PP + (M0*DELTA_T)*(AA+GG) + ruidos

DPP(i,j,k)=(PPNEW-phase_field_old(i,j,k))

phase_field_new(i,j,k)=dmax1(0.d0,dmin1(1.d0,PPNEW))
