subroutine concentracao_soluto

use variaveis
implicit none
CHARACTER (LEN=5)  :: tecform


! Loop-invariant reciprocals -- scalars, captured by value at each target launch.
InvDX=1.D0/(DX*DX)
InvDY=1.D0/(DY*DY)
InvDZ=1.D0/(DZ*DZ)


if (use_gpu) then

  ! -------------------------------------------------------------------------
  ! GPU path. phase_field_old, conc_soluto_old/new resident on the device.
  ! Kernel 1 -- device copy new<=old (carries halo cells the interior sweep skips).
  ! -------------------------------------------------------------------------
  !$omp target teams distribute parallel do collapse(3) private(i,j,k)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        conc_soluto_new(i,j,k) = conc_soluto_old(i,j,k)
      end do
    end do
  end do
  !$omp end target teams distribute parallel do

  ! -------------------------------------------------------------------------
  ! Kernel 2 -- solute-diffusion stencil. Read-only: phase_field_old,
  ! conc_soluto_old. Written at owned (i,j,k) only: conc_soluto_new. Neighbour
  ! reads are all from conc_soluto_old / phase_field_old (double buffer) -> no
  ! loop-carried dependence, every (i,j,k) independent. PRIVATE list identical to
  ! the CPU directive.
  ! -------------------------------------------------------------------------
  !$omp target teams distribute parallel do collapse(3) &
  !$omp private(i,j,k) &
  !$omp private(PP,PN,PS,PW,PE,PT,PB) &
  !$omp private(CPP,CPS,CPN,CPW,CPE,CPT,CPB) &
  !$omp private(h,h_E,h_W,h_N,h_S,h_T,h_B) &
  !$omp private(CSP,CSN,CSS,CSE,CSW,CST,CSB) &
  !$omp private(CLP,CLN,CLS,CLE,CLW,CLT,CLB) &
  !$omp private(D_soluto_C) &
  !$omp private(GAMA_P,GAMA_S,GAMA_N,GAMA_W,GAMA_E,GAMA_T,GAMA_B) &
  !$omp private(GAMAs,GAMAn,GAMAw,GAMAe,GAMAt,GAMAb) &
  !$omp private(ran,ruidos,c_new)
  do K=2, nz-1
  do J=2, ny-1
   do I=2, nx-1
     include 'conc_body.inc'
   end do
  end do
  end do
  !$omp end target teams distribute parallel do

  ! Kernel 3 -- device copy old<=new (interior), matching the original serial
  ! "conc_soluto_old = conc_soluto_new" so the boundary decoration below can read
  ! the just-updated interior from conc_soluto_old.
  !$omp target teams distribute parallel do collapse(3) private(i,j,k)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        conc_soluto_old(i,j,k) = conc_soluto_new(i,j,k)
      end do
    end do
  end do
  !$omp end target teams distribute parallel do

  ! Kernel 4 -- boundary decoration on the device (single launch, same array
  ! sections as the CPU path).
  !$omp target
  conc_Soluto_new(1:NX,1,1:NZ)     = Conc_Soluto_Old(1:NX,2,1:NZ)
  conc_Soluto_new(1:NX,ny,1:NZ)    = Conc_Soluto_Old(1:NX,ny-1,1:NZ)
  Conc_Soluto_new(1,1:NY,1:NZ)     = Conc_Soluto_Old(2,1:NY,1:NZ)
  Conc_Soluto_new(nx,1:NY,1:NZ)    = Conc_Soluto_Old(nx-1,1:NY,1:NZ)
  Conc_Soluto_new(1:NX,1:NY,1)     = Conc_Soluto_Old(1:NX,1:NY,2)
  Conc_Soluto_new(1:NX,1:NY,nz)    = Conc_Soluto_Old(1:NX,1:NY,nz-1)
  !$omp end target

else

! ===========================================================================
! CPU path -- UNCHANGED original OpenMP implementation.
! ===========================================================================
conc_soluto_new = conc_soluto_old

! Already double-buffered in the original (writes only conc_soluto_new during the sweep,
! swaps into conc_soluto_old only after the loop). Same DEFAULT(NONE)/full-PRIVATE-list
! rationale as the other subroutines: every scalar below is a `variaveis` module-global.
!$omp parallel default(none) &
!$omp shared(phase_field_old,conc_soluto_old,conc_soluto_new) &
!$omp shared(nx,ny,nz,DELTA_T,InvDX,InvDY,InvDZ,Ke,AAA,D_liq,D_sol,it) &
!$omp private(i,j,k) &
!$omp private(PP,PN,PS,PW,PE,PT,PB) &
!$omp private(CPP,CPS,CPN,CPW,CPE,CPT,CPB) &
!$omp private(h,h_E,h_W,h_N,h_S,h_T,h_B) &
!$omp private(CSP,CSN,CSS,CSE,CSW,CST,CSB) &
!$omp private(CLP,CLN,CLS,CLE,CLW,CLT,CLB) &
!$omp private(D_soluto_C) &
!$omp private(GAMA_P,GAMA_S,GAMA_N,GAMA_W,GAMA_E,GAMA_T,GAMA_B) &
!$omp private(GAMAs,GAMAn,GAMAw,GAMAe,GAMAt,GAMAb) &
!$omp private(ran,ruidos,c_new)
!$omp do collapse(2) schedule(static)
do K=2, nz-1
do J=2, ny-1
 do I=2, nx-1
   include 'conc_body.inc'
end do
	end do
	end do
!$omp end do
!$omp end parallel

	conc_soluto_old = conc_soluto_new

 !=========  Decora��o do Contorno Inclu�do no C�lculo ==============================
	 conc_Soluto_new(1:NX,1,1:NZ)     = Conc_Soluto_Old(1:NX,2,1:NZ)
	 conc_Soluto_new(1:NX,ny,1:NZ)    = Conc_Soluto_Old(1:NX,ny-1,1:NZ)
	 Conc_Soluto_new(1,1:NY,1:NZ)     = Conc_Soluto_Old(2,1:NY,1:NZ)
     Conc_Soluto_new(nx,1:NY,1:NZ)    = Conc_Soluto_Old(nx-1,1:NY,1:NZ)
     Conc_Soluto_new(1:NX,1:NY,1)     = Conc_Soluto_Old(1:NX,1:NY,2)
     Conc_Soluto_new(1:NX,1:NY,nz)    = Conc_Soluto_Old(1:NX,1:NY,nz-1)

end if
!===================================================================================


end subroutine concentracao_soluto
