!ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
   Subroutine Calc_Temperature()
!ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
   USE variaveis
   Implicit None
   Integer Iter
   Double Precision,Allocatable :: AW(:,:), AE(:,:), AS(:,:), AN(:,:), ap0(:,:), ap(:,:), sp_(:,:), su(:,:)
   Double Precision :: RES
   Double Precision :: Tol = 1.d+05

!==================================================================================================
 

 allocate (AW(NX, NY), AE(NX,NY), AN(NX,NY), AS(NX,NY), AP(NX,NY), sp_(NX,NY), su(NX,NY), ap0(nx,ny))

DO iter=1,IT_MAX
 Call Coeficientes ()
 Call TempBound()
 Call Source()
 Call Solver()
 IF(DABS(RES)<=Tol) EXIT
ENDDO

Return 

Contains 

Subroutine Coeficientes()
Implicit None 
double precision :: D, hp, PP

!-- Initializing Coefficients:

 do i = 1, nx
    do j = 1, ny

       PP = Phase_Field_new(I,J)
       hp = PP
!       hp = PP**3.d0*(6.d0*PP**2.d0 - 15.0*PP - 10.d0)
       D  = D_s * hp + (1.d0 - hp) * D_l
       AW(i,j)= D*DY/DX
       AE(i,j)= D*DY/DX
       AS(i,j)= D*DX/DY
       AN(i,j)= D*DX/DY
    end do
end do

End subroutine


Subroutine TempBound()
Implicit None 
Integer :: I,J
!---- Boundary Conditions: 
   

	DO I=1, nx 
	   Do j=1, ny 
!	   T_(i,1)  = T_(i,2)    ; AS(i,2)   = 0.d0
       T_(i,1)  = T0 - 20.0  ; AS(i,2)   = 0.d0
       T_(i,ny) = T_(i,ny-1) ; AN(i,ny-1)= 0.d0
       T_(1,j)  = T_(2,j)    ; AW(2,j)   = 0.d0
	   T_(nx,j) = T_(nx-1,j) ; AE(nx-1,j)= 0.d0
       Enddo
	ENDDO

	
	

End subroutine TempBound

! ==================================================
! Montagem dos termos fontes e do coeficiente AP
! ==================================================
Subroutine Source()
Implicit None 
Integer :: I,J
Double Precision::Hlinha, PP, hp, cp

    DO i=2,NX-1
	   Do j=2,NY-1
			ap0(i,j)=DX*DY/DELTA_t
			sp_(i,j)=-ap0(i,j)
			PP=Phase_Field_new(I,J)
            Hlinha  = (3.0*(PP**2.0))*(10.0-15.0*PP+6.0*(PP**2.0))+(PP**3.0)*(-15.0+12.0*PP)
            hp = PP**3.d0*(6.d0*PP**2.d0 - 15.0*PP - 10.d0)
            cp  = cp_s + (1.d0 - hp) * cp_l
            su(i,j)= ap0(i,j)*TOLD(i,j)+(DELTA_H/cp)*(DX*DY/DELTA_t)*Hlinha*(phase_field_new(I,J)-phase_field_old(I,J))
			ap(i,j)=AW(i,j)+AE(i,j)+AS(i,j)+AN(i,j)-sp_(i,j)
	   End do
	End do
	return

End subroutine Source 



Subroutine Solver()
Implicit None 

Integer :: I,J
Double Precision :: Bet,temp, TNew_
Double Precision,Allocatable ::  qt(:),pt(:)
    



!   i-direction 
    ALLOCATE(qt(1:NX),pt(1:NX))
	RES=0.d0
    DO j=2,NY-1
      pt(1)=0.d0
      qt(1)=T_(1,j)

      DO i=2,NX-1
        bet=ap(i,j)-pt(i-1)*AW(i,j)
        IF(DABS(bet) <= small)bet=small
        pt(i)=AE(i,j)/bet
        temp=AS(i,j)*T_(i,j-1)+AN(i,j)*T_(i,j+1)+su(i,j)
        qt(i)=(temp+AW(i,j)*qt(i-1))/bet
      END DO
      DO i=NX-1,2,-1
	    Tnew_=T_(i+1,j)*pt(i)+qt(i)
	    RES=DMAX1(RES,DABS((Tnew_-T_(i,j))/T_(i,j)))
        T_(i,j)=Tnew_
      END DO
    END DO
    DEALLOCATE(qt,Pt)

!--- J - direction:
    ALLOCATE(qt(1:NY),pt(1:NY))
    RES=0.d0
    DO i=2,NX-1
      pt(1)=0.d0
      qt(1)=T_(i,1)
      DO j=2,NY-1
        bet=ap(i,j)-pt(j-1)*AS(i,j)
        IF(DABS(bet) <= small)bet=small
        pt(j)=AN(i,j)/bet
        temp=AW(i,j)*T_(i-1,j)+AE(i,j)*T_(i+1,j)+su(i,j)
        qt(j)=(temp+AS(i,j)*qt(j-1))/bet
      END DO
      DO j=NY-1,2,-1
	    Tnew_=T_(i,j+1)*pt(j)+qt(j)
	    RES=DMAX1(RES,DABS((Tnew_-T_(i,j))/T_(i,j)))
        T_(i,j)=Tnew_
      END DO
    END DO
    DEALLOCATE(qt,Pt)

!--- I - direction, reversed:
    ALLOCATE(qt(1:NX),pt(1:NX))
    RES=0.d0
    DO j=NY-1,2,-1
      pt(1)=0.d0
      qt(1)=T_(1,j)
      DO i=2,NX-1
        bet=ap(i,j)-pt(i-1)*AW(i,j)
        IF(DABS(bet) <= small)bet=small
        pt(i)=AE(i,j)/bet
        temp=AS(i,j)*T_(i,j-1)+AN(i,j)*T_(i,j+1)+su(i,j)
        qt(i)=(temp+AW(i,j)*qt(i-1))/bet
      END DO
      DO i=NX-1,2,-1
	    Tnew_=T_(i+1,j)*pt(i)+qt(i)
	    RES=DMAX1(RES,DABS((Tnew_-T_(i,j))/T_(i,j)))
        T_(i,j)=Tnew_
      END DO
    END DO

    DEALLOCATE(qt,Pt)

 !--- J - direction, reversed :
    ALLOCATE(qt(1:NY),pt(1:NY)) 
    RES=0.d0
    DO i=NX-1,2,-1
      pt(1)=0.d0
      qt(1)=T_(i,1)
      DO j=2,NY-1
        bet=ap(i,j)-pt(j-1)*AS(i,j)
        IF(DABS(bet) <= small)bet=small
        pt(j)=AN(i,j)/bet
        temp=AW(i,j)*T_(i-1,j)+AE(i,j)*T_(i+1,j)+su(i,j)
        qt(j)=(temp+AS(i,j)*qt(j-1))/bet
      END DO
      DO j=NY-1,2,-1
	    Tnew_=T_(i,j+1)*pt(j)+qt(j)
	    RES=DMAX1(RES,DABS((Tnew_-T_(i,j))/T_(i,j)))
        T_(i,j)=Tnew_
      END DO
    END DO
    DEALLOCATE(qt,Pt)

End subroutine Solver


End Subroutine  Calc_Temperature
!ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ


