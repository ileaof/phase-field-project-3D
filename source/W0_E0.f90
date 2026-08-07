Subroutine W0_E0
Use Variaveis
Implicit none
Double precision::Epslon_Try, Epslon





DO EPSLON_Try=0,1,0.0000001


        W0 = ((1.5556349186104d0/Lambda)*Epslon_Try)*((1.5556349186104d0/Lambda)*Epslon_Try)


        E0 = ajuste*Sigma*4.24264068711929d0/sqrt(W0)


        write(*,*)"E0-EPSLON_Try==>", abs(E0-EPSLON_Try)


        IF (ABS(E0-EPSLON_Try)<1D-05) THEN
                go to 10
        End if



END DO


10 WRITE(*,*) "W0==>",W0,"     ","E0==>",E0,'     ',"Dt==>",delta_t



End Subroutine W0_E0
