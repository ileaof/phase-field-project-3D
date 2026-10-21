subroutine entrada

use variaveis

implicit none
double precision :: Teta0_inic
integer :: ios_tp
! Condicao inicial: tipo_inicial = 0 -> frente plana em j<=4 (original);
!                   tipo_inicial = 1 -> semente esferica no centro do dominio
!                   tipo_inicial = 2 -> 1/8 de semente no canto
!                   tipo_inicial = 3 -> frente plana na base (k<=4), crescimento em +z
integer :: tipo_inicial, ic, jc, kc
double precision :: R_semente, dist2

!open (unit=1,file='/Users/ileao/Documents/Phase_field/Al-Cu/source/entrada.dat')
open (unit=1,file='entrada.dat')

read (1,*)
read (1,*)
read (1,*)
read (1,*) Tm
read (1,*) T0
read (1,*) delta_H
read (1,*) E0
read (1,*) w0
read (1,*) M0
read (1,*) D_s
read (1,*) D_l
read (1,*) Cp_s
read (1,*) Cp_l
read (1,*) Nx
read (1,*) Ny
read (1,*) Nz
read (1,*) dy
read (1,*) dx
read (1,*) dz
read (1,*) small
read (1,*) ANI
read (1,*) it_max
read (1,*) AAA
read (1,*) Teta0_inic
read (1,*) DELTA_EPSLON
read (1,*) DELTA_T
read (1,*) D_sol
read (1,*) D_liq
read (1,*) Ke
read (1,*) Me
read (1,*) Vm
read (1,*) R
read (1,*) C0
read (1,*) Thin_Interface
read (1,*) Sigma
read (1,*) Beta
read (1,*) ajuste
read (1,*) stoptime
read (1,*) WriteFileD_timeCounter
! TP: temperatura isotermica usada nas equacoes de fase e soluto. Linha opcional
! para manter compatibilidade com entrada.dat antigos: se faltar, usa 919 K
! (o valor que antes ficava fixo em principal.f90).
read (1,*,iostat=ios_tp) TP
if (ios_tp /= 0) then
    TP = 919.0d0
    write(*,"('AVISO: linha TP ausente no entrada.dat -- usando TP = 919.0 K')")
end if
! Linha opcional: tipo de condicao inicial e raio da semente (em celulas).
! Se faltar, mantem a frente plana original.
read (1,*,iostat=ios_tp) tipo_inicial, R_semente
if (ios_tp /= 0) then
    tipo_inicial = 0
    R_semente    = 0.d0
end if
! Linha opcional: taxa de resfriamento [K/s] e varreduras da fase por passo,
! como no codigo 2D (TP = T0 - taxa*t; laco Do N=1,2). Se faltar: TP fixo, 1 varredura.
read (1,*,iostat=ios_tp) taxa_resfr, n_varreduras
if (ios_tp /= 0) then
    taxa_resfr   = 0.d0
    n_varreduras = 1
end if
n_varreduras = max(1, n_varreduras)
if (taxa_resfr > 0.d0) TP = T0


close (unit=1)

write(*,"(/,'### Dados de Entrada ###',/)")
write(*,"('Tm = ',F14.2,' K')")Tm 
write(*,"('T0 = ',F14.2,' K')")T0
write(*,"('delta_H = ',E14.6,' J/m^3')")delta_H 
write(*,"('E0 = ',E14.6,' -')")E0 
write(*,"('w0 = ',E14.6,' -')")w0 
write(*,"('M0 = ',E14.6,' -')")M0
write(*,"('D_s = ',E14.6,' m^2/s')")D_s 
write(*,"('D_l = ',E14.6,' m^2/s')")D_l 
write(*,"('Cp_s = ',E14.6,' J/kg.K')")Cp_s 
write(*,"('Cp_l = ',E14.6,' J/kg.K')")Cp_l 
write(*,"('Nx = ',I5,' -')")Nx
write(*,"('Ny = ',I5,' -')")Ny
write(*,"('Nz = ',I5,' -')")Nz
write(*,"('dy = ',E14.6,' m')")dy
write(*,"('dx = ',E14.6,' m')")dx
write(*,"('dz = ',E14.6,' m')")dz
write(*,"('small = ',E14.6,' -')")small 
write(*,"('ANI = ',E14.6,' -')") ANI
write(*,"('it_max = ',I10,' -')") it_max
write(*,"('ANI = ',E14.6,' -')") ANI
write(*,"('AAA = ',F14.4,' -')") AAA
write(*,"('Theta0_inic = ',F14.4,' -')")Teta0_inic
write(*,"('DELTA_EPSLON = ',F14.4,' -')")DELTA_EPSLON
write(*,"('DELTA_t = ',E14.4,' s')")DELTA_t
write(*,"('D_sol = ',E14.4,' m^2/s')")D_sol
write(*,"('D_liq = ',E14.4,' m^2/s')")D_liq
write(*,"('Ke = ',F14.4,' -')")Ke
write(*,"('Me = ',F14.4,' K/mol')")Me
write(*,"('Mv = ',E14.4,' m^3/mol')")Vm
write(*,"('R = ',F14.5,' J/mol K')")R
write(*,"('C0 = ',F14.5,' molar fraction')")C0
write(*,"('Thin_Interface = ',F14.3,' molar fraction')")Thin_Interface
write(*,"('Sigma = ',F14.3,' J/m^2')")Sigma
write(*,"('ajuste',F14.3,' -')")ajuste
write(*,"('stoptime',F14.3,' s')")stoptime
write(*,"('WriteFileD_timeCounter',I8,' -')")WriteFileD_timeCounter
write(*,"('TP = ',F14.2,' K  (Cl_eq = ',F10.6,', TL(C0) = ',F10.2,' K)')") TP, (Tm-TP)/Me, Tm-Me*C0
if (taxa_resfr > 0.d0) then
    write(*,"('Resfriamento continuo: TP = T0 - ',F8.3,' K/s * t  (T0 = ',F8.2,' K)')") taxa_resfr, T0
end if
write(*,"('Varreduras da fase por passo = ',I0)") n_varreduras
if (tipo_inicial == 3) then
    write(*,"('Condicao inicial: frente plana na base (k<=4), crescimento em +z')")
else if (tipo_inicial == 2) then
    write(*,"('Condicao inicial: 1/8 de semente no canto (3,3,3), raio = ',F6.1,' celulas')") R_semente
else if (tipo_inicial == 1) then
    write(*,"('Condicao inicial: semente esferica central, raio = ',F6.1,' celulas')") R_semente
else
    write(*,"('Condicao inicial: frente plana (j<=4)')")
end if

write(*,"(/,'### Fim Leitura de Entrada de Dados ###',/)")


!Ke = 2*Ke  xxxxxxx     BOA FIGURA PARA O TRABALHO  xxxxxxxxxxxxxxxxxxxxxxxxxx
! Ke = 0.001*Ke
! D_liq=2D0*D_liq

!Me = 1.2d0*Me


!============CALCULO DA ESPESSURA DA MALHA========================================================
   XL1 = DX*NX ; YL1= DY*NY ; ZL1 = DZ*NZ
!============CALCULO DA ESPESSURA DA INTERFACE====================================================
    LAMBDA = Thin_Interface*dx
! Explicit-diffusion stability limit is dt<=dx^2/(4D) in 2D (2 second-derivative terms)
! but dt<=dx^2/(6D) in 3D (3 terms); the 2D code's divisor of 5 (80% of the 2D limit of 4)
! is scaled to 7.5 (80% of the 3D limit of 6) to keep the same stability margin in 3D.
    !DELTA_T = (DX*DX)/(7.5d0*D_LIQ)
    !DELTA_T = (DX*DX)/(10.d0*D_LIQ)
    !DELTA_T = (DX*DX)/(20.d0*D_LIQ) !Data Sim 00
    !DELTA_T = (DX*DX)/(30.d0*D_LIQ) !Data Sim 01
    !DELTA_T = (DX*DX)/(40.d0*D_LIQ)  !Data Sim 02
    !DELTA_T = (DX*DX)/(50.d0*D_LIQ)  !Data Sim 03
    !DELTA_T = (DX*DX)/(60.d0*D_LIQ)  !Data Sim 03
    !DELTA_T = (DX*DX)/(12.d0*D_LIQ)  ! metade do limite 3D (1/6); 5x mais rapido que /60
    DELTA_T = (DX*DX)/(7.5d0*D_LIQ)  ! equivalente 3D do /5 do codigo 2D (80% do limite)
   



!========CALCULO DO DT PARA CONDI��O DE CONVERGENCIA==============================================
  !DELTA_T=0.3*((XL1/NX)*(XL1/NX))/D   !CRIT�RIO DE CONVERG�NCIA    
!============CALCULO DA ESPESSURA DA MALHA========================================================




!call random(ran_0)
!TETA_0 (crystal-axis rotation) dropped in 3D: cubic anisotropy is expressed directly via
!grad(phi) components, no single angle. Teta0_inic is still read above for entrada.dat
!backward-compatibility but is otherwise unused (it was always 0 in practice, and the only
!code that varied it per-grain was already commented out).




!=======================Dimensionalizando os vetores================================

allocate(phase_field_old(1:Nx, 1:Ny, 1:Nz),phase_field_new(1:Nx, 1:Ny, 1:Nz))
allocate( T_(1:Nx,1:Ny,1:Nz))
allocate(xc(1:Nx, 1:Ny, 1:Nz),yc(1:Nx, 1:Ny, 1:Nz),zc(1:Nx, 1:Ny, 1:Nz),DPP(1:Nx, 1:Ny, 1:Nz),TOLD(1:Nx, 1:Ny, 1:Nz),TNEW(1:Nx, 1:Ny, 1:Nz))
allocate(conc_soluto_new(1:Nx, 1:Ny, 1:Nz), conc_soluto_old(1:Nx, 1:Ny, 1:Nz))



!============================Zerando os vetores (first-touch paralelo)===============
! Initialized via the SAME parallel loop pattern (K,J,I nesting, collapse(2), static
! schedule) as the main solver loops, so each thread "first-touches" (and becomes the
! NUMA-local owner of) exactly the pages it will repeatedly read/write during time
! stepping, instead of a single serial thread mapping every page to one NUMA node.
! DPP is also zeroed here: temperatura.f90 reads DPP(I,J,K) before equacao_campo_fase.f90
! ever writes it on the very first iteration (principal.f90 calls temperatura before
! equacao_campo_fase), so it was previously read uninitialized on IT=1 -- closed as a
! safe, behavior-preserving-in-intent fix (zero is what an unwritten DPP should mean:
! no phase-field change yet).
ic = Nx/2 ; jc = Ny/2 ; kc = Nz/2
! tipo_inicial = 2: semente no canto (1/8 de semente). As paredes de fluxo zero
! funcionam como planos de simetria, dobrando o alcance de cada braco.
if (tipo_inicial == 2) then
    ic = 3 ; jc = 3 ; kc = 3
end if

!$omp parallel default(none) &
!$omp shared(phase_field_old,phase_field_new,TNEW,TOLD,conc_soluto_new,conc_soluto_old,DPP) &
!$omp shared(T0,C0,Ke,Nx,Ny,Nz,tipo_inicial,R_semente,ic,jc,kc) private(i,j,k,dist2)
!$omp do collapse(2) schedule(static)
do k=1,Nz
  do j=1,Ny
    do i=1,Nx
      if (tipo_inicial == 1 .or. tipo_inicial == 2) then
        dist2 = dble(i-ic)**2 + dble(j-jc)**2 + dble(k-kc)**2
      else
        dist2 = huge(1.d0)
      end if
      if ((tipo_inicial == 0 .and. j<=4) .or. (tipo_inicial == 3 .and. k<=4) &
          .or. dist2 <= R_semente**2) then
        phase_field_old(i,j,k) = 1.d0
        conc_soluto_old(i,j,k) = Ke*C0
      else
        phase_field_old(i,j,k) = 0.d0
        conc_soluto_old(i,j,k) = C0
      end if
      phase_field_new(i,j,k) = 0.d0
      TNEW(i,j,k) = T0
      TOLD(i,j,k) = T0
      conc_soluto_new(i,j,k) = C0
      DPP(i,j,k) = 0.d0
    end do
  end do
end do
!$omp end do
!$omp end parallel

!phase_field_old(0.80*nx,0.80*ny,0.80*nz)=+1.d0
!conc_soluto_old(0.80*nx,0.80*ny,0.80*nz)=Ke*C0


!phase_field_old(0.40*nx,0.40*ny,0.40*nz)=+1.d0
!conc_soluto_old(0.40*nx,0.40*ny,0.40*nz)=Ke*C0


TX = 10
!==================================================================================


!==================================================================================


end subroutine entrada

