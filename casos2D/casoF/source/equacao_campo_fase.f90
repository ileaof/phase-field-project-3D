subroutine equacao_campo_fase
use variaveis
implicit none



	 phase_field_new=phase_field_old


     !av = 50+0.2d-02*delta_T*IT


Do N=1,2

IF (N==1) THEN
D1 = 0.60*nx
!==============nucleo setado=========================================================================
!!phase_field_old(0.60*nx,0.60*ny)  =    +1.d0
!!conc_soluto_old(0.60*nx,0.60*ny)  =    Ke*C0
!call random(ran_0)
!!TETA_0(D1-TX:D1+TX,D1-TX:D1+TX,1)               =    0.D0   !ran_0  
END IF

IF (N==2) THEN
D2=0.30*nx
!==============nucleo setado=========================================================================
!!phase_field_old(0.30*nx,0.30*ny)  =    +1.d0
!!conc_soluto_old(0.30*nx,0.30*ny)  =    Ke*C0
!call random(ran_0)
!!TETA_0(D2-TX:D2+TX,D2-TX:D2+TX,2)               =    0.785398163d0  ! ran_0  
END IF


Do J=3, Ny-2
 Do I=3, nx-2



!====================================IDENTIFICAÇÃO DOS PONTOS========================================

PP = phase_field_old(i,j)

PN= phase_field_old(i,j+1)

PS= phase_field_old(i,j-1)

PE= phase_field_old(i+1,j)

PW= phase_field_old(i-1,j)

PNW= phase_field_old(i-1,j+1)

PSW= phase_field_old(i-1,j-1)

PNE= phase_field_old(i+1,j+1)

PSE= phase_field_old(i+1,j-1)

PPNEW= phase_field_new(i,j) 

!TP = Told(i,j)
      
!==============================================PONTO E===============================================

small =1.d-23


PHIy_E          =  (PNE - PSE)*0.5d0/dy
PHIx_E          =  (phase_field_old(i+2,j)-PP)*0.5/dx
TETA_E          =   DATAN(PHIy_E/(PHIx_E + small))
EPSLON_E_LINHA  =  -E0*DELTA_EPSLON*ANI*DSIN(ANI*(TETA_E - TETA_0(I,J,N))) 
EPSLON_E        =   E0*(1+DELTA_EPSLON*DCOS(ANI*(TETA_E - TETA_0(I,J,N))))



!=============================================PONTO W================================================

PHIy_W          =  (PNW-PSW)*0.5d0/dy
PHIx_W          =  (PP-phase_field_old (i-2,j))*0.5d0/dx
TETA_W          =   DATAN(PHIy_W/(PHIx_W+small))
EPSLON_W_LINHA  =  -E0*DELTA_EPSLON*ANI*DSIN(ANI*(TETA_W - TETA_0(I,J,N))) 
EPSLON_W        =   E0*(1+DELTA_EPSLON*DCOS(ANI*(TETA_W - TETA_0(I,J,N))))



!============================================PONTO N=================================================

PHIy_N          =  (phase_field_old(i,j+2)-PP)*0.5d0/dy
PHIx_N          =  (PNE-PNW)*0.5d0/dx
TETA_N          =   DATAN(PHIy_N/(PHIx_N + small))
EPSLON_N_LINHA  =  -E0*DELTA_EPSLON*ANI*DSIN(ANI*(TETA_N - TETA_0(I,J,N))) 
EPSLON_N        =   E0*(1+DELTA_EPSLON*DCOS(ANI*(TETA_N - TETA_0(I,J,N))))



!============================================PONTO S=================================================

PHIy_S          =  (PP-phase_field_old (i,j-2))*0.5d0/dy
PHIx_S          =  (PSE-PSW)*0.5d0/dx
TETA_S          =   DATAN(PHIy_S/(PHIx_S + small))
EPSLON_S_LINHA  =  -E0*DELTA_EPSLON*ANI*DSIN(ANI*(TETA_S - TETA_0(I,J,N))) 
EPSLON_S        =   E0*(1+DELTA_EPSLON*DCOS(ANI*(TETA_S - TETA_0(I,J,N))))



!=========================================EQUAÇÃO DE FASE============================================


AA = (0.25/(dx*dy))*(PNE*(EPSLON_N*EPSLON_N_LINHA - EPSLON_E*EPSLON_E_LINHA)+&
                     PSE*(EPSLON_E*EPSLON_E_LINHA - EPSLON_S*EPSLON_S_LINHA)+&
					 PNW*(EPSLON_W*EPSLON_W_LINHA - EPSLON_N*EPSLON_N_LINHA)+&
					 PSW*(EPSLON_S*EPSLON_S_LINHA - EPSLON_W*EPSLON_W_LINHA))
					 
BB = PN*((0.25d0/(dx*dy))*(EPSLON_W*EPSLON_W_LINHA-EPSLON_E*EPSLON_E_LINHA)+((EPSLON_N**2.d0)/(dy**2.d0)))				  


CC = PS*((0.25d0/(dx*dy))*(EPSLON_E*EPSLON_E_LINHA-EPSLON_W*EPSLON_W_LINHA)+((EPSLON_S**2.d0)/(dy**2.d0)))


DD = PE*((0.25d0/(dx*dy))*(EPSLON_N*EPSLON_N_LINHA-EPSLON_S*EPSLON_S_LINHA)+((EPSLON_E**2.d0)/(dx**2.d0)))


EE = PW*((0.25d0/(dx*dy))*(EPSLON_S*EPSLON_S_LINHA-EPSLON_N*EPSLON_N_LINHA)+((EPSLON_W**2.d0)/(dx**2.d0)))


FF = PP*(-(1/(dx**2.d0))*((EPSLON_E**2.d0)+(EPSLON_W**2.d0))-(1/(dy**2.d0))*((EPSLON_N**2.d0)+(EPSLON_S**2.d0))+(1.d0/(M0*DELTA_T)))



!====================================================================================================

g_linha = 2.d0*PP*(1.d0-PP)*(1.d0-PP) - 2.d0*PP*PP*(1.d0-PP)                     !PP*(2-(6*PP)+(4*(PP**2))) 

h_linha = 30.d0*PP*PP*(1.d0-PP)*(1.d0-PP) !  !(3.D0*PP*PP)*(6.D0*PP*PP-15.D0*PP+10.D0) + (PP*PP*PP)*(12.D0*PP-15.D0)    !(PP**2)*(30-(60*PP)+(30*(PP**2))) 

CPP  = conc_soluto_old(i,j)

h   = PP*PP*PP*(10.d0-15.d0*PP+6.d0*PP*PP) !   !PP*PP*PP*(6.d0*PP*PP-15.d0*PP+10.d0)

CLP  = CPP/(Ke*h+(1.d0-h))

CSP  = Ke*CPP/(Ke*h+(1.d0-h))   

Cs_eq = ((Tm-TP)*Ke)/ Me

Cl_eq = (Tm-TP)/ Me


	  


GG  = (-w0*g_linha )+ (R*TP/Vm)*h_linha*Log(abs(((1.d0-CS_eq)*(1.d0-CLP))/((1.d0-CL_eq)*(1.d0-CSP))))   


!========================================Geracao de ruidos===========================================

 call random(ran)

 ruidos= (2.d0*ran-1.d0)*AAA*16.d0*(PP**2)*(1-PP)**2.d0



PPNEW = (M0*DELTA_T)*(AA+BB+CC+DD+EE+FF+GG) + ruidos

DPP(i,j)=(PPNEW-phase_field_old(i,j))


IF (.false. .and. ABS(phase_field_old(0.30*nx+TX,0.30*ny+TX))>1.0D-06)THEN
TX = TX+20
END IF


phase_field_new(i,j)=dmax1(0.d0,dmin1(1.d0,PPNEW))

phase_field_old(i,j)=phase_field_new(i,j)

end do

end do
 end do

!=========  Decoração do Contorno Incluído no Cálculo ==============================
	 Phase_Field_new(1:NX,1)   = Phase_Field_new(1:NX,2)
	 Phase_Field_new(1:NX,2)   = Phase_Field_new(1:NX,3)

	 Phase_Field_new(1:NX,ny)  = Phase_Field_new(1:NX,ny-1)
	 Phase_Field_new(1:NX,ny-1)  = Phase_Field_new(1:NX,ny-2)

	 Phase_Field_new(1,1:NY)   = Phase_Field_new(2,1:NY)
	 Phase_Field_new(2,1:NY)   = Phase_Field_new(3,1:NY)

     Phase_Field_new(nx,1:NY)  = Phase_Field_new(nx-1,1:NY)
     Phase_Field_new(nx-1,1:NY)  = Phase_Field_new(nx-2,1:NY)

!==================================================================================================

write(*,4635) phase_field_old(nx/2,10),tp, phase_field_old(nx/2,100)                             !0.30*ny+TX),0.30*nx-TX,0.30*nx+TX,0.30*ny-TX,0.30*ny+TX
                             !0.30*ny+TX),0.30*nx-TX,0.30*nx+TX,0.30*ny-TX,0.30*ny+TX
!TETA_0(nx/2,ny/2,1),TETA_0(nx/2,ny/2,2),TX,0.70*nx+TX           !ABS(DPP(0.70*nx,0.75*ny)),av  !M0*E0*E0                        !,TETA_1,TETA_2, TETA_3

! ,T0-Tp,DELTA_T*IT    !       GG   
!DELTA_T*IT 
4635  FORMAT(9(1X,e14.4))


end subroutine  equacao_campo_fase