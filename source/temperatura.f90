subroutine temperatura

use variaveis

implicit none

real*8 cp

cp = cp_s

! ===========================================================================
! GPU path: TOLD/TNEW/phase_field_old/DPP are already resident on the device
! (mapped once by principal.f90's `target enter data`), so nothing crosses the
! PCIe/NVLink bus here. Two target kernels: (1) carry the full field forward
! (TNEW=TOLD, needed for the halo cells the interior loop skips); (2) the
! stencil update. collapse(3) flattens the whole K,J,I nest so the runtime can
! spread all interior points across teams/threads (maximum occupancy) with the
! fastest index i mapped to consecutive threads (coalesced global-memory loads).
! ===========================================================================
if (use_gpu) then

  !$omp target teams distribute parallel do collapse(3) &
  !$omp private(i,j,k)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        TNEW(i,j,k) = TOLD(i,j,k)
      end do
    end do
  end do
  !$omp end target teams distribute parallel do

  !$omp target teams distribute parallel do collapse(3) &
  !$omp private(i,j,k,TE,TW,TN,TS,TT,TB,PP,h_linha,FONTE,TNEW1)
  DO K=2,NZ-1
     DO J=2,NY-1
     DO I=2,NX-1
       include 'temp_body.inc'
     END DO
     END DO
  END DO
  !$omp end target teams distribute parallel do

else

! ===========================================================================
! CPU path -- UNCHANGED from the original OpenMP implementation.
! Told is read-only inside this loop (double-buffered Jacobi update -- the former
! in-place "TOLD(I,J,K)=TNEW(I,J,K)" write has been removed, matching the same fix
! applied to equacao_campo_fase.f90; see accompanying analysis). IMPORTANT:
! principal.f90's end-of-iteration "TNEW = TOLD" was relying on that in-place write
! to have already made TOLD equal to TNEW by this point (making that line a harmless
! no-op); it has been changed to "TOLD = TNEW" there to correctly advance the field
! in time now that this subroutine no longer does it internally.
TNEW=TOLD

!$omp parallel default(none) &
!$omp shared(Told,TNEW,phase_field_old,DPP,nx,ny,nz,dx,dy,dz,D,DELTA_T,delta_H,cp,TP) &
!$omp private(i,j,k,TE,TW,TN,TS,TT,TB,PP,h_linha,FONTE,TNEW1)
!$omp do collapse(2) schedule(static)
DO K=2,NZ-1
     DO J=2,NY-1
     DO I=2,NX-1
       include 'temp_body.inc'
     END DO
     END DO
  END DO
!$omp end do
!$omp end parallel

end if


!WRITE(*,*)phase_field_old(nx*0.9,ny*0.9), TOLD(nx*0.9,ny*0.9), it




END SUBROUTINE temperatura
