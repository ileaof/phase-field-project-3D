subroutine temperatura

use variaveis

implicit none

real*8 cp

cp = cp_s

TNEW=TOLD



10  DO J=2,NY-1
     DO I=2,NX-1

!TP= Told(i,j)

TE= Told(i+1,j)

TW= Told(i-1,j)

TN= Told(i,j+1)

TS= Told(i,j-1)


      
	   PP = phase_field_old(I,J)


!=======================================TERMO FONTE======================================


  h_linha = 30.d0*PP*PP*(1.d0-PP)*(1.d0-PP) !        ! (ileao mesmo do fi) estava este --> !(PP**2)*(30-(60*PP)+(30*(PP**2)))  
 FONTE   = (delta_H/Cp)*h_linha*DPP(I,J)             ! (DELTAH/cp)*(DX*DY/Dtime)*Hlinha*(Phase_Field(I,J)-Phase_Field_Old(I,J))


!========================================================================================

!TNEW1 = TOLD(I,J)+((D*DELTA_T)/(DX*DX))*(TOLD(I+1,J)-4*TOLD(I,J)+TOLD(I-1,J)+TOLD(I,J+1)+TOLD(I,J-1))  + FONTE

TNEW1= D*DELTA_T*(((TE-2.d0*TP+TW)/dx**2.d0)+((TN-2.d0*TP+TS)/dy**2.d0))+FONTE+TP


     
	 TNEW(I,J) = TNEW1
	 TOLD(I,J)=TNEW(I,J)


 END DO
  END DO


!WRITE(*,*)phase_field_old(nx*0.9,ny*0.9), TOLD(nx*0.9,ny*0.9), it    




END SUBROUTINE temperatura