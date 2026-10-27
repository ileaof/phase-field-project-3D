!ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
   Subroutine Saida !(JTime)
!ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
   USE Variaveis
   Implicit None 
   Integer :: JTime 
   Double Precision :: xp, yp, C1, C2
   CHARACTER (LEN=5)  :: tecform
   CHARACTER (LEN=100) :: OUTPUTFILENAME
   CHARACTER (LEN=52) :: Path_to_output
!==================================================================================================
  
!  Select Case (jTime) 
!	Case (1)   ! Arquivo para campos iniciais 
!		  OUTPUTFILENAME="saida_inicial.dat"
!	Case default
!		  OUTPUTFILENAME="saida.dat"
! End select 

  write (text,'(I6)') IT
  write (text2,'(f14.7)') time
  text3 = adjustl('saida_'//trim(text)//'_time_'//trim(text2)//'s'//'.dat')
!  Path_to_output = "/Users/ileao/Documents/Phase_field/Al-Cu/output/"
  Path_to_output = ".\output\"
  OUTPUTFILENAME = adjustl(Path_to_output//trim(text3))
!  OPEN(UNIT=2,FILE= "saida2017.dat" ,STATUS="unknown")     !OUTPUTFILENAME,STATUS="unknown")
  OPEN(UNIT=2, FILE= OUTPUTFILENAME, STATUS="unknown")     !OUTPUTFILENAME,STATUS="unknown")
  tecform="POINT"
  WRITE(2,*)'TITLE="Phase Field"'
  WRITE(2,4535)
  WRITE(2,*)"ZONE"
  WRITE(2,10001)NX,NY,tecform
  IF(tecform == "POINT") THEN
    
    DO j=1,NY
      DO i=1,NX
          xp=MAX(i-1,0)*DX
          yp=MAX(j-1,0)*DY
          WRITE(2,4635)Xp,  Yp, phase_field_new(I,J),TP,Conc_Soluto_New(I,J) !TNEW(I,J)
      END DO
    END DO
  END IF
  WRITE(2,*)"       "
  WRITE(2,*)"       "
  CLOSE(UNIT=2)
  4535   FORMAT('variables= "X","Y","PhaseField","Temp","Conc_Solute"')
  4635   FORMAT(5(1X,e14.7))
 
10001  FORMAT( 1X,'T=" " ,I=',i4,',J=',i4,',F=',a5)
!WRITE(*,*)
!WRITE(*,*)" +-----------------------------------------+"
!WRITE(*,*)" |   Tecplot files has been generated......|"
!WRITE(*,*)" +-----------------------------------------+"
!WRITE(*,*)


!==================================================================================================
   End Subroutine  Saida
!ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
