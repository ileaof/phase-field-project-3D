subroutine diagnosticos(L2_phi, phi_max, phi_min, solute_total, free_energy)
! Verification/reduction diagnostics for comparing serial vs. OpenMP-parallel vs. GPU runs.
! None of the main solver loops contain a reduction (every grid point's update is
! independent), so this is the one place that legitimately needs
! reduction(+:...)/reduction(max:...)/reduction(min:...).
use variaveis
implicit none
double precision, intent(out) :: L2_phi, phi_max, phi_min, solute_total, free_energy
double precision :: phi_p, grad2, w_term

L2_phi = 0.d0
phi_max = -huge(1.d0)
phi_min = huge(1.d0)
solute_total = 0.d0
free_energy = 0.d0

if (use_gpu) then

  ! GPU path. phase_field_new / conc_soluto_new are resident on the device, so only
  ! the five scalar accumulators cross the bus (map(tofrom:...) seeds them with the
  ! host initial values above and copies the final reduced results back). The GPU
  ! runtime performs an efficient tree reduction across teams and threads.
  !$omp target teams distribute parallel do collapse(3) &
  !$omp map(tofrom:L2_phi,solute_total,free_energy,phi_max,phi_min) &
  !$omp private(i,j,k,phi_p,grad2,w_term) &
  !$omp reduction(+:L2_phi,solute_total,free_energy) &
  !$omp reduction(max:phi_max) reduction(min:phi_min)
  do k=1,Nz
    do j=1,Ny
      do i=1,Nx
        include 'diag_body.inc'
      end do
    end do
  end do
  !$omp end target teams distribute parallel do

else

  ! CPU path -- UNCHANGED original OpenMP reduction.
  !$omp parallel do collapse(2) schedule(static) default(none) &
  !$omp shared(phase_field_new,conc_soluto_new,Nx,Ny,Nz,dx,dy,dz,E0,w0) &
  !$omp private(i,j,k,phi_p,grad2,w_term) &
  !$omp reduction(+:L2_phi,solute_total,free_energy) &
  !$omp reduction(max:phi_max) reduction(min:phi_min)
  do k=1,Nz
    do j=1,Ny
      do i=1,Nx
        include 'diag_body.inc'
      end do
    end do
  end do
  !$omp end parallel do

end if

L2_phi = sqrt(L2_phi)

end subroutine diagnosticos
