subroutine entrada

use variaveis

implicit none
double precision :: Teta0_inic

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
read (1,*) dy
read (1,*) dx
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


close (unit=1)

write(*,"(/,'### Dados de Entrada ###',/)")
write(*,"('Tm = ',F12.2,' K')")Tm 
write(*,"('T0 = ',F12.2,' K')")T0
write(*,"('delta_H = ',E12.6,' J/m^3')")delta_H 
write(*,"('E0 = ',E12.6,' -')")E0 
write(*,"('w0 = ',E12.6,' -')")w0 
write(*,"('M0 = ',E12.6,' -')")M0
write(*,"('D_s = ',E12.6,' m^2/s')")D_s 
write(*,"('D_l = ',E12.6,' m^2/s')")D_l 
write(*,"('Cp_s = ',E12.6,' J/kg.K')")D_s 
write(*,"('Cp_l = ',E12.6,' J/kg.K')")D_l 
write(*,"('Nx = ',I5,' -')")Nx 
write(*,"('Ny = ',I5,' -')")Ny 
write(*,"('dy = ',E12.6,' m')")dy 
write(*,"('dx = ',E12.6,' m')")dx 
write(*,"('small = ',E12.6,' -')")small 
write(*,"('ANI = ',E12.6,' -')") ANI
write(*,"('it_max = ',I10,' -')") it_max
write(*,"('ANI = ',E12.6,' -')") ANI
write(*,"('AAA = ',F12.4,' -')") AAA
write(*,"('Theta0_inic = ',F12.4,' -')")Teta0_inic
write(*,"('DELTA_EPSLON = ',F12.4,' -')")DELTA_EPSLON
write(*,"('DELTA_t = ',E12.4,' s')")DELTA_t
write(*,"('D_sol = ',E12.4,' m^2/s')")D_sol
write(*,"('D_liq = ',E12.4,' m^2/s')")D_liq
write(*,"('Ke = ',F12.4,' -')")Ke
write(*,"('Me = ',F12.4,' K/mol')")Me
write(*,"('Mv = ',E12.4,' m^3/mol')")Vm
write(*,"('R = ',F12.5,' J/mol K')")R
write(*,"('C0 = ',F12.5,' molar fraction')")C0
write(*,"('Thin_Interface = ',F12.3,' molar fraction')")Thin_Interface
write(*,"('Sigma = ',F12.3,' J/m^2')")Sigma
write(*,"('ajuste',F12.3,' -')")ajuste
write(*,"('stoptime',F12.3,' s')")stoptime
write(*,"('WriteFileD_timeCounter',I4,' -')")WriteFileD_timeCounter

write(*,"(/,'### Fim Leitura de Entrada de Dados ###',/)")


!Ke = 2*Ke  xxxxxxx     BOA FIGURA PARA O TRABALHO  xxxxxxxxxxxxxxxxxxxxxxxxxx
! Ke = 0.001*Ke
! D_liq=2D0*D_liq

!Me = 1.2d0*Me


!============CALCULO DA ESPESSURA DA MALHA========================================================
   XL1 = DX*NX ; YL1= DY*NY
!============CALCULO DA ESPESSURA DA INTERFACE====================================================
    LAMBDA = Thin_Interface*dx
    DELTA_T = (DX*DX)/(5.d0*D_LIQ)
    !DELTA_T = (DX*DX)/(10*D_LIQ)




!========CALCULO DO DT PARA CONDIÇÃO DE CONVERGENCIA==============================================
  !DELTA_T=0.3*((XL1/NX)*(XL1/NX))/D   !CRITÉRIO DE CONVERGÊNCIA    
!============CALCULO DA ESPESSURA DA MALHA========================================================




!call random(ran_0)
!TETA_0          =   0.d0   !ran_0  
!!call random(ran_1)
!!TETA_1          =   0.785398163d0 !ran_1  
!!call random(ran_2)
!!TETA_2          =   ran_2  
!!call random(ran_3)
!!TETA_3          =   ran_3  




!=======================Dimensionalizando os vetores================================

allocate(phase_field_old(1:Nx, 1:Ny),phase_field_new(1:Nx, 1:Ny), TETA_0(1:NX,1:NY,1:2))
allocate( T_(1:Nx,1:Ny))
allocate(xc(1:Nx, 1:Ny),yc(1:Nx, 1:Ny),DPP(1:Nx, 1:Ny),TOLD(1:Nx, 1:Ny),TNEW(1:Nx, 1:Ny))
allocate(conc_soluto_new(1:Nx, 1:Ny), conc_soluto_old(1:Nx, 1:Ny))
teta_0=Teta0_inic
!allocate(phase_field_old(1:Nx, 1:1000),phase_field_new(1:Nx, 1:1000))
!allocate(xc(1:Nx, 1:1000),yc(1:Nx, 1:1000),DPP(1:Nx, 1:1000),TOLD(1:Nx, 1:1000),TNEW(1:Nx, 1:1000))
!allocate(conc_soluto_new(1:Nx, 1:1000), conc_soluto_old(1:Nx, 1:1000))



!============================Zerando os vetores=====================================


phase_field_old= +0.d0
phase_field_new= +0.d0

TNEW     = T0
TOld     = TNEW

conc_soluto_new = C0
conc_soluto_old = C0

!=============================NUCLEO e CONCENTRAÇÃO INICIAL==========================
phase_field_old(1:NX,1:4)=+1.d0
conc_soluto_old(1:NX,1:4)=Ke*C0

!phase_field_old(0.80*nx,0.80*ny)=+1.d0
!conc_soluto_old(0.80*nx,0.80*ny)=Ke*C0


!phase_field_old(0.40*nx,0.40*ny)=+1.d0
!conc_soluto_old(0.40*nx,0.40*ny)=Ke*C0


TX = 10
!==================================================================================


!==================================================================================


end subroutine entrada

