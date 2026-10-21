Program  principal

use variaveis

implicit none

call entrada
call W0_E0


time=0.d0

Do it=1, it_max


            time=delta_T+time

            call temperatura
            !!call Calc_Temperature  ! Estava comentada
            !!TNEW = T_

            TP = T0 - 19.84d0*DELTA_T*IT          !0.005*DELTA_T*IT   !83.3333D0

            !========================================================================================
            !DIF_CONC=ABS(C_NEW-C0)
            !nyf=ny
            !IF(DIF_CONC>=1.D-02) THEN
            !NYF=NY+1
            !NY=NYF
            !END IF

            !WRITE(*,*)NYF,DIF_CONC

            !========================================================================================


            !IF (DELTA_T*it>=16.2D-02) THEN     !  7.920D-02) THEN
            !WRITE(*,*)'FIM DOS CALCULOS'
            !STOP
            !END IF


            call equacao_campo_fase

            call concentracao_soluto
            
            !Test of local equilibrium statement
            DO I = 2, NX-1
                DO J = 2, NY-1
                   if (conc_soluto_new(I,J).GT. C0) then
                             if (phase_field_new(I,J) .LT. 0.9d0) then
                                   TNEW(I,J) = Tm - Me*conc_soluto_new(I,J) 
                             endif
                   endif
                END DO
            END DO

            !IF(1.0d+13>=IT>=1) THEN
            call auxiliar
            !END IF


            if(E_R >= WriteFileD_timeCounter) then

            E_R = 0

                      WRITE(*,2000)
                      WRITE(*,2010)
                      WRITE(*,2000)
                      WRITE(*,2000)
                      WRITE(*,2020) IT
                      WRITE(*,2030) time
                      WRITE(*,2040) stoptime
                      WRITE(*,2050) WriteFileD_timeCounter*Delta_T
 !                     WRITE(*,*)'ESCREVENDO OS RESULTADOS ','  ','TETA_0====',ran_0,'  ','TETA_1====',ran_1 
 !                     WRITE(*,*)'ESCREVENDO OS RESULTADOS ','  ','TETA_2====',ran_2,'  ','TETA_3====',ran_3 
                      WRITE(*,2000)
            call saida
            end if
            e_r = e_r + 1



            phase_field_old = phase_field_new
            TNEW   = TOLD
            conc_soluto_old = conc_soluto_new
            

            !if(IT>= 8000) then
            !go to 100
            !end if

            if ((time .GE. (stoptime)).AND.(stoptime .NE. 0.d0)) then
                      write(*,1000) stoptime
                      WRITE(*,2000)
                      WRITE(*,2010)
                      WRITE(*,2000)
                      WRITE(*,2000)
                      WRITE(*,2020) IT
                      WRITE(*,2030) time
                      WRITE(*,2040) stoptime
                      WRITE(*,2050) WriteFileD_timeCounter*Delta_T
 !                     WRITE(*,*)'ESCREVENDO OS RESULTADOS ','  ','TETA_0====',ran_0,'  ','TETA_1====',ran_1 
 !                     WRITE(*,*)'ESCREVENDO OS RESULTADOS ','  ','TETA_2====',ran_2,'  ','TETA_3====',ran_3 
                      WRITE(*,2000)
                      if(E_R < WriteFileD_timeCounter ) call saida
                      exit
             endif
                
            


end do

!! 100 call saida   ! última saida repetida by ileao ! must be included a Delta_T on the stoptime eq. above

!============== liberando os vetores da memória alocada =====================

deallocate(phase_field_old,phase_field_new, TETA_0)
deallocate(xc,yc,DPP,TOLD,TNEW)
deallocate(conc_soluto_new, conc_soluto_old)
deallocate( T_ )

!=============================================================================
1000  format (//,5X,'Program terminated by stop time condition = ',f14.7,' s',/)
2000  format ('=====================================================')
2010  format ('=============  ESCREVENDO RESULTADOS  ===============')
2015  format ('Parada por criterio de tempo total < tempo de escrita')
2020  format ('IT==> ', I8,' iteracoes')
2030  format ('Tempo total simulado ==> ',d12.4,' s')
2040  format ('Tempo de parada (stoptime) ==> ', d12.4,' s')
2050  format ('Intervalo de escrita em arquivo  ==> ', d12.4,' s')

end program  principal