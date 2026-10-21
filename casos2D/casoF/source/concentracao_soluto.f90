subroutine concentracao_soluto

use variaveis
implicit none
CHARACTER (LEN=5)  :: tecform


InvDX=1.D0/(DX*DX)
InvDY=1.D0/(DY*DY)

conc_soluto_new = conc_soluto_old

do J=2, ny-1
 do I=2, nx-1

!===========================IDENTIFICAÇÃO DAS CONCENTRAÇÕES=============================

 PP = phase_field_old(i,j)

 PN = phase_field_old(i,j+1)

 PS = phase_field_old(i,j-1)

 PW = phase_field_old(i-1,j)

 PE = phase_field_old(i+1,j)


CPP = conc_soluto_old(i,j)

CPS = conc_soluto_old(i,j-1)

CPN = conc_soluto_old(i,j+1)

CPW = conc_soluto_old(i-1,j)

CPE = conc_soluto_old(i+1,j)


  h = (PP**3)*(10.d0-(15.d0*PP)+(6.d0*(PP**2)))   !(PP**3)*(10.d0-(15.d0*PP)+(6.d0*(PP**2)))

h_E = (PE**3)*(10.d0-(15.d0*PE)+(6.d0*(PE**2)))   !(PE**3)*(10.d0-(15.d0*PE)+(6.d0*(PE**2)))

h_W = (PW**3)*(10.d0-(15.d0*PW)+(6.d0*(PW**2)))   !(PW**3)*(10.d0-(15.d0*PW)+(6.d0*(PW**2)))

h_N = (PN**3)*(10.d0-(15.d0*PN)+(6.d0*(PN**2)))   !(PN**3)*(10.d0-(15.d0*PN)+(6.d0*(PN**2)))

h_S = (PS**3)*(10.d0-(15.d0*PS)+(6.d0*(PS**2)))   !(PS**3)*(10.d0-(15.d0*PS)+(6.d0*(PS**2)))


CSP = (Ke*CPP)/((Ke*h)+(1.d0-h))

CSN = (Ke*CPN)/((Ke*h_N)+(1.d0-h_N))

CSS = (Ke*CPS)/((Ke*h_S)+(1.d0-h_S))

CSE = (Ke*CPE)/((Ke*h_E)+(1.d0-h_E))

CSW = (Ke*CPW)/((Ke*h_W)+(1.d0-h_W))


CLP = CPP/(Ke*h+(1.d0-h))

CLN = CPN/(Ke*h_N+(1.d0-h_N))

CLS = CPS/(Ke*h_S+(1.d0-h_S))

CLE = CPE/(Ke*h_E+(1.d0-h_E))

CLW = CPW/(Ke*h_W+(1.d0-h_W))

!====================================================================================================
!IF ( (I<=0.5*NX).AND.(J<=0.5*NY) ) THEN
!Teta_0 = 0.785398163D0
!phase_field_old(0.25*NX,0.25*NY)=+1.d0
!conc_soluto_old(0.25*NX,0.25*NY)=Ke*C0
!END IF

!IF ( (I>0.5*NX).AND.(J<=0.5*NY) ) THEN
!Teta_0 = 0.0D0
!phase_field_old(0.75*NX,0.25*NY)=+1.d0
!conc_soluto_old(0.75*NX,0.25*NY)=Ke*C0

!END IF

!IF ( (I<=0.5*NX).AND.(J>0.5*NY) ) THEN
!Teta_0 = 0.785398163D0
!phase_field_old(0.25*NX,0.75*NY)=+1.d0
!conc_soluto_old(0.25*NX,0.75*NY)=Ke*C0

!END IF

!IF ( (I>0.5*NX).AND.(J>0.5*NY) ) THEN
!Teta_0 = 0.0D0
!phase_field_old(0.75*NX,0.75*NY)=+1.d0
!conc_soluto_old(0.75*NX,0.75*NY)=Ke*C0

!END IF
!==============================================PONTO E===============================================
!IF(c_new>=0.0196) THEN
!phase_field_old(i,j)=0.D0
!END IF
!====================CONDIÇÃO DE DIFUSIVIDADE NO SOLIDO E NO LIQUIDO====================

!IF (PP<=0.1) THEN
!	D_soluto_C = D_liq
!END IF

!IF (PP>=0.9) THEN
!	D_soluto_C = D_sol
!END IF

!IF (0.1<PP<0.9) THEN
!	D_soluto_C = D_liq
!END IF


If (PP<=0.95d0)then      !  0.87d0

	D_soluto_C = D_liq

Else

	D_soluto_C = D_sol

End if


!===================IDENTIFICAÇÃO DOS GAMA'S USADOS NA EQUAÇÃO DE CONCENTRAÇÃO===========


GAMA_P = D_soluto_C *((1.d0-h)*CLP*(1.d0-CLP)+h*CSP*(1.d0-CSP))

GAMA_S = D_soluto_C *((1.d0-h_S)*CLS*(1.d0-CLS)+h_S*CSS*(1.d0-CSS))

GAMA_N = D_soluto_C *((1.d0-h_N)*CLN*(1.d0-CLN)+h_N*CSN*(1.d0-CSN))

GAMA_W = D_soluto_C *((1.d0-h_W)*CLW*(1.d0-CLW)+h_W*CSW*(1.d0-CSW))

GAMA_E = D_soluto_C *((1.d0-h_E)*CLE*(1.d0-CLE)+h_E*CSE*(1.d0-CSE))


GAMAs = 2.d0/ ((1.d0/GAMA_P)+(1.d0/GAMA_S))

GAMAn = 2.d0/ ((1.d0/GAMA_N)+(1.d0/GAMA_P))

GAMAw = 2.d0/ ((1.d0/GAMA_P)+(1.d0/GAMA_W))

GAMAe = 2.d0/ ((1.d0/GAMA_E)+(1.d0/GAMA_P))

!========================================================================================
call random(ran)

!IF (PP<0.001) THEN
ruidos= (1)*ran*AAA*16*(phase_field_old(i,j)**2)*(1.d0-phase_field_old(i,j))**2
!ELSE
!RUIDOS=0.D0
!END IF


	c_new = (GAMAe*InvDX*(log(CLE/(1.d0-CLE)) - log(CLP/(1.d0-CLP)))- &
	         GAMAw*InvDX*(log(CLP/(1.d0-CLP)) - log(CLW/(1.d0-CLW)))+ &
			 GAMAn*InvDY*(log(CLN/(1.d0-CLN)) - log(CLP/(1.d0-CLP)))- &
			 GAMAs*InvDY*(log(CLP/(1.d0-CLP)) - log(CLS/(1.d0-CLS))))*DELTA_T + CPP ! + DELTA_T*RUIDOS

!========================================================================================


	conc_soluto_new(i,j) = c_new


end do
	end do

	conc_soluto_old = conc_soluto_new


 !=========  Decoração do Contorno Incluído no Cálculo ==============================
	 conc_Soluto_new(1:NX,1)   = Conc_Soluto_Old(1:NX,2)
	 conc_Soluto_new(1:NX,ny)  = Conc_Soluto_Old(1:NX,ny-1)
	 Conc_Soluto_new(1,1:NY)   = Conc_Soluto_Old(2,1:NY)
     Conc_Soluto_new(nx,1:NY)  = Conc_Soluto_Old(nx-1,1:NY)

!===================================================================================
!if((iT*DELTA_T)>=5.5D-05) then  !estava como comentario
!stop
!write(*,*)"fim", "m0==>",m0
!end if
!===================================================================================





!conc_soluto_new(nx/2, 10+ny/2), T0-Tp   !conc_soluto_new(nx/2, 10),,E0,W0 


end subroutine concentracao_soluto	
				  