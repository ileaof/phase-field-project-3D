! ============================================================================
! conc_body.inc -- per-cell body of the solute-diffusion update
! ----------------------------------------------------------------------------
! Single source of truth, INCLUDEd verbatim from both the CPU worksharing loop
! and the GPU offload loop in concentracao_soluto.f90. Every name written below
! is a loop index or a module-global scalar listed PRIVATE on both enclosing
! directives, so each (i,j,k) iteration owns its own copy -- race-free on CPU
! threads and GPU teams/threads alike.
! ============================================================================

 PP = phase_field_old(i,j,k)

 PN = phase_field_old(i,j+1,k)

 PS = phase_field_old(i,j-1,k)

 PW = phase_field_old(i-1,j,k)

 PE = phase_field_old(i+1,j,k)

 PT = phase_field_old(i,j,k+1)

 PB = phase_field_old(i,j,k-1)


CPP = conc_soluto_old(i,j,k)

CPS = conc_soluto_old(i,j-1,k)

CPN = conc_soluto_old(i,j+1,k)

CPW = conc_soluto_old(i-1,j,k)

CPE = conc_soluto_old(i+1,j,k)

CPT = conc_soluto_old(i,j,k+1)

CPB = conc_soluto_old(i,j,k-1)


! h(phi) quintico, o mesmo da equacao de fase (phase_body.inc). Anterior: h = PP*PP*(3-2*PP)
  h = (PP**3)*(10.d0-(15.d0*PP)+(6.d0*(PP**2)))

h_E = (PE**3)*(10.d0-(15.d0*PE)+(6.d0*(PE**2)))

h_W = (PW**3)*(10.d0-(15.d0*PW)+(6.d0*(PW**2)))

h_N = (PN**3)*(10.d0-(15.d0*PN)+(6.d0*(PN**2)))

h_S = (PS**3)*(10.d0-(15.d0*PS)+(6.d0*(PS**2)))

h_T = (PT**3)*(10.d0-(15.d0*PT)+(6.d0*(PT**2)))

h_B = (PB**3)*(10.d0-(15.d0*PB)+(6.d0*(PB**2)))


CSP = (Ke*CPP)/((Ke*h)+(1.d0-h))

CSN = (Ke*CPN)/((Ke*h_N)+(1.d0-h_N))

CSS = (Ke*CPS)/((Ke*h_S)+(1.d0-h_S))

CSE = (Ke*CPE)/((Ke*h_E)+(1.d0-h_E))

CSW = (Ke*CPW)/((Ke*h_W)+(1.d0-h_W))

CST = (Ke*CPT)/((Ke*h_T)+(1.d0-h_T))

CSB = (Ke*CPB)/((Ke*h_B)+(1.d0-h_B))


CLP = CPP/(Ke*h+(1.d0-h))

CLN = CPN/(Ke*h_N+(1.d0-h_N))

CLS = CPS/(Ke*h_S+(1.d0-h_S))

CLE = CPE/(Ke*h_E+(1.d0-h_E))

CLW = CPW/(Ke*h_W+(1.d0-h_W))

CLT = CPT/(Ke*h_T+(1.d0-h_T))

CLB = CPB/(Ke*h_B+(1.d0-h_B))

!====================CONDI��O DE DIFUSIVIDADE NO SOLIDO E NO LIQUIDO====================

If (PP<=0.95d0)then      !  0.87d0

	D_soluto_C = D_liq

Else

	D_soluto_C = D_sol

End if


!===================IDENTIFICA��O DOS GAMA'S USADOS NA EQUA��O DE CONCENTRA��O===========


GAMA_P = D_soluto_C *((1.d0-h)*CLP*(1.d0-CLP)+h*CSP*(1.d0-CSP))

GAMA_S = D_soluto_C *((1.d0-h_S)*CLS*(1.d0-CLS)+h_S*CSS*(1.d0-CSS))

GAMA_N = D_soluto_C *((1.d0-h_N)*CLN*(1.d0-CLN)+h_N*CSN*(1.d0-CSN))

GAMA_W = D_soluto_C *((1.d0-h_W)*CLW*(1.d0-CLW)+h_W*CSW*(1.d0-CSW))

GAMA_E = D_soluto_C *((1.d0-h_E)*CLE*(1.d0-CLE)+h_E*CSE*(1.d0-CSE))

GAMA_T = D_soluto_C *((1.d0-h_T)*CLT*(1.d0-CLT)+h_T*CST*(1.d0-CST))

GAMA_B = D_soluto_C *((1.d0-h_B)*CLB*(1.d0-CLB)+h_B*CSB*(1.d0-CSB))


GAMAs = 2.d0/ ((1.d0/GAMA_P)+(1.d0/GAMA_S))

GAMAn = 2.d0/ ((1.d0/GAMA_N)+(1.d0/GAMA_P))

GAMAw = 2.d0/ ((1.d0/GAMA_P)+(1.d0/GAMA_W))

GAMAe = 2.d0/ ((1.d0/GAMA_E)+(1.d0/GAMA_P))

GAMAt = 2.d0/ ((1.d0/GAMA_T)+(1.d0/GAMA_P))

GAMAb = 2.d0/ ((1.d0/GAMA_P)+(1.d0/GAMA_B))

!========================================================================================
ran = pseudo_random(i,j,k,it,2)

ruidos= (1)*ran*AAA*16*(phase_field_old(i,j,k)**2)*(1.d0-phase_field_old(i,j,k))**2


	c_new = (GAMAe*InvDX*(log(CLE/(1.d0-CLE)) - log(CLP/(1.d0-CLP)))- &
	         GAMAw*InvDX*(log(CLP/(1.d0-CLP)) - log(CLW/(1.d0-CLW)))+ &
			 GAMAn*InvDY*(log(CLN/(1.d0-CLN)) - log(CLP/(1.d0-CLP)))- &
			 GAMAs*InvDY*(log(CLP/(1.d0-CLP)) - log(CLS/(1.d0-CLS)))+ &
			 GAMAt*InvDZ*(log(CLT/(1.d0-CLT)) - log(CLP/(1.d0-CLP)))- &
			 GAMAb*InvDZ*(log(CLP/(1.d0-CLP)) - log(CLB/(1.d0-CLB))))*DELTA_T + CPP ! + DELTA_T*RUIDOS

!========================================================================================


	conc_soluto_new(i,j,k) = c_new
