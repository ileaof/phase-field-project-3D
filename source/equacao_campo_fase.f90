subroutine equacao_campo_fase
use variaveis
implicit none



! Note: the 2D code's outer "Do N=1,2" loop existed only to index a per-grain crystal
! orientation TETA_0(:,:,N). In 3D the cubic anisotropy is expressed directly from
! grad(phi) (no single orientation angle), TETA_0 no longer exists, and the grain-specific
! nucleation code that would have distinguished N=1 from N=2 was already commented out in
! the original (both passes were physically identical) -- so that loop is dropped and this
! runs a single pass.

D1 = 0.60*nx
D2 = 0.30*nx

! Cs_eq/Cl_eq depend only on TP,Tm,Ke,Me -- NOT on (i,j,k) -- loop-invariant. Computed
! once here on the HOST. Because they are plain scalars, each target launch below captures
! them by value (implicit firstprivate), so the freshly computed values reach the device
! automatically; auxiliar.f90 (host) also reads them afterward to compute M0.
Cs_eq = ((Tm-TP)*Ke)/ Me
Cl_eq = (Tm-TP)/ Me


if (use_gpu) then

  ! -------------------------------------------------------------------------
  ! GPU path. phase_field_old/new, conc_soluto_old, DPP are resident on the
  ! device (mapped once in principal.f90). No host<->device transfer occurs in
  ! the hot loop.
  !
  ! Kernel 1 -- device copy new<=old. Fills the halo/boundary cells that the
  ! interior sweep (3..N-2) never writes, exactly as the CPU "phase_field_new=
  ! phase_field_old" whole-array assignment does.
  ! -------------------------------------------------------------------------
  !$omp target teams distribute parallel do collapse(3) private(i,j,k)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        phase_field_new(i,j,k) = phase_field_old(i,j,k)
      end do
    end do
  end do
  !$omp end target teams distribute parallel do

  ! -------------------------------------------------------------------------
  ! Kernel 2 -- the anisotropic phase-field stencil. Read-only: phase_field_old,
  ! conc_soluto_old. Written at the owned (i,j,k) only: phase_field_new, DPP.
  ! Every neighbour access is a pure READ of phase_field_old (Jacobi / double
  ! buffer), so there is no loop-carried dependence -> every (i,j,k) is
  ! independent and safe to run on distinct teams/threads. collapse(3) maximises
  ! exposed parallelism and keeps i (unit stride) on consecutive threads for
  ! coalesced loads. The long PRIVATE list is IDENTICAL to the CPU directive so
  ! the two paths compute the same values; each iteration gets its own scratch.
  ! -------------------------------------------------------------------------
  !$omp target teams distribute parallel do collapse(3) &
  !$omp private(i,j,k) &
  !$omp private(PP,PN,PS,PE,PW,PT,PB) &
  !$omp private(PNE,PSE,PNW,PSW,PNT,PNB,PST,PSB,PET,PEB,PWT,PWB) &
  !$omp private(PPNEW,small,K_aniso) &
  !$omp private(phix_E,phiy_E,phiz_E,Q_E,S4_E,eps_E,exy_E,ezx_E) &
  !$omp private(phix_W,phiy_W,phiz_W,Q_W,S4_W,eps_W,exy_W,ezx_W) &
  !$omp private(phix_N,phiy_N,phiz_N,Q_N,S4_N,eps_N,exy_N,eyz_N) &
  !$omp private(phix_S,phiy_S,phiz_S,Q_S,S4_S,eps_S,exy_S,eyz_S) &
  !$omp private(phix_T,phiy_T,phiz_T,Q_T,S4_T,eps_T,eyz_T,ezx_T) &
  !$omp private(phix_B,phiy_B,phiz_B,Q_B,S4_B,eps_B,eyz_B,ezx_B) &
  !$omp private(AA,Cross_xy,Cross_yz,Cross_zx) &
  !$omp private(g_linha,h_linha,CPP,h,CLP,CSP,GG) &
  !$omp private(ran,ruidos)
  Do K=3, Nz-2
  Do J=3, Ny-2
   Do I=3, nx-2
     include 'phase_body.inc'
   end do
  end do
  end do
  !$omp end target teams distribute parallel do

  ! -------------------------------------------------------------------------
  ! Boundary decoration on the device. Same array-section assignments as the CPU
  ! path, wrapped in ONE target region (single kernel launch, executed on the
  ! device so phase_field_new never has to come back to the host). These planes
  ! are O(N^2) -- negligible next to the O(N^3) stencil.
  ! -------------------------------------------------------------------------
  !$omp target
  Phase_Field_new(1:NX,1,1:NZ)     = Phase_Field_new(1:NX,2,1:NZ)
  Phase_Field_new(1:NX,2,1:NZ)     = Phase_Field_new(1:NX,3,1:NZ)
  Phase_Field_new(1:NX,ny,1:NZ)    = Phase_Field_new(1:NX,ny-1,1:NZ)
  Phase_Field_new(1:NX,ny-1,1:NZ)  = Phase_Field_new(1:NX,ny-2,1:NZ)
  Phase_Field_new(1,1:NY,1:NZ)     = Phase_Field_new(2,1:NY,1:NZ)
  Phase_Field_new(2,1:NY,1:NZ)     = Phase_Field_new(3,1:NY,1:NZ)
  Phase_Field_new(nx,1:NY,1:NZ)    = Phase_Field_new(nx-1,1:NY,1:NZ)
  Phase_Field_new(nx-1,1:NY,1:NZ)  = Phase_Field_new(nx-2,1:NY,1:NZ)
  Phase_Field_new(1:NX,1:NY,1)     = Phase_Field_new(1:NX,1:NY,2)
  Phase_Field_new(1:NX,1:NY,2)     = Phase_Field_new(1:NX,1:NY,3)
  Phase_Field_new(1:NX,1:NY,nz)    = Phase_Field_new(1:NX,1:NY,nz-1)
  Phase_Field_new(1:NX,1:NY,nz-1)  = Phase_Field_new(1:NX,1:NY,nz-2)
  !$omp end target

  ! Pull back ONLY the two cells the progress print below reads (a few bytes),
  ! keeping the full field resident on the device.
  !$omp target update from( phase_field_new(nx/2,10,max(nz/2,1)), &
  !$omp                     phase_field_new(nx/2,min(100,ny-1),max(nz/2,1)) )

else

! ===========================================================================
! CPU path -- UNCHANGED original OpenMP implementation.
! ===========================================================================
	 phase_field_new=phase_field_old

! Parallelized over (K,J) via COLLAPSE(2), leaving I as the innermost serial/SIMD loop
! (Fortran is column-major: I is the fastest-varying array index, so keeping it innermost
! and un-split gives each thread large contiguous, cache-line-aligned write chunks in
! phase_field_new -- no false sharing at thread boundaries). Every scalar touched inside
! the loop body is module-global (declared in `variaveis`, not subroutine-local), so ALL
! of them must be listed PRIVATE below or every thread would clobber the same variable;
! DEFAULT(NONE) forces a compile error instead of a silent race if one is missed.
! phase_field_old is READ-ONLY inside this loop (double-buffered Jacobi update).
!$omp parallel default(none) &
!$omp shared(phase_field_old,phase_field_new,conc_soluto_old,DPP) &
!$omp shared(nx,ny,nz,dx,dy,dz,E0,DELTA_EPSLON,DELTA_T,M0,w0,Tm,Me,Vm,R,Ke,AAA,TP,it) &
!$omp shared(Cs_eq,Cl_eq) &
!$omp private(i,j,k) &
!$omp private(PP,PN,PS,PE,PW,PT,PB) &
!$omp private(PNE,PSE,PNW,PSW,PNT,PNB,PST,PSB,PET,PEB,PWT,PWB) &
!$omp private(PPNEW,small,K_aniso) &
!$omp private(phix_E,phiy_E,phiz_E,Q_E,S4_E,eps_E,exy_E,ezx_E) &
!$omp private(phix_W,phiy_W,phiz_W,Q_W,S4_W,eps_W,exy_W,ezx_W) &
!$omp private(phix_N,phiy_N,phiz_N,Q_N,S4_N,eps_N,exy_N,eyz_N) &
!$omp private(phix_S,phiy_S,phiz_S,Q_S,S4_S,eps_S,exy_S,eyz_S) &
!$omp private(phix_T,phiy_T,phiz_T,Q_T,S4_T,eps_T,eyz_T,ezx_T) &
!$omp private(phix_B,phiy_B,phiz_B,Q_B,S4_B,eps_B,eyz_B,ezx_B) &
!$omp private(AA,Cross_xy,Cross_yz,Cross_zx) &
!$omp private(g_linha,h_linha,CPP,h,CLP,CSP,GG) &
!$omp private(ran,ruidos)
!$omp do collapse(2) schedule(static)
Do K=3, Nz-2
Do J=3, Ny-2
 Do I=3, nx-2
   include 'phase_body.inc'
end do

end do
end do
!$omp end do
!$omp end parallel

!=========  Decora��o do Contorno Inclu�do no C�lculo ==============================
	 Phase_Field_new(1:NX,1,1:NZ)     = Phase_Field_new(1:NX,2,1:NZ)
	 Phase_Field_new(1:NX,2,1:NZ)     = Phase_Field_new(1:NX,3,1:NZ)

	 Phase_Field_new(1:NX,ny,1:NZ)    = Phase_Field_new(1:NX,ny-1,1:NZ)
	 Phase_Field_new(1:NX,ny-1,1:NZ)  = Phase_Field_new(1:NX,ny-2,1:NZ)

	 Phase_Field_new(1,1:NY,1:NZ)     = Phase_Field_new(2,1:NY,1:NZ)
	 Phase_Field_new(2,1:NY,1:NZ)     = Phase_Field_new(3,1:NY,1:NZ)

     Phase_Field_new(nx,1:NY,1:NZ)    = Phase_Field_new(nx-1,1:NY,1:NZ)
     Phase_Field_new(nx-1,1:NY,1:NZ)  = Phase_Field_new(nx-2,1:NY,1:NZ)

     Phase_Field_new(1:NX,1:NY,1)     = Phase_Field_new(1:NX,1:NY,2)
     Phase_Field_new(1:NX,1:NY,2)     = Phase_Field_new(1:NX,1:NY,3)

     Phase_Field_new(1:NX,1:NY,nz)    = Phase_Field_new(1:NX,1:NY,nz-1)
     Phase_Field_new(1:NX,1:NY,nz-1)  = Phase_Field_new(1:NX,1:NY,nz-2)

end if
!==================================================================================================

! Progress print. On the GPU path the two sampled cells were just pulled back with a tiny
! `target update from` above; tp is a host scalar. (print-only, not a physics change).
write(*,4635) phase_field_new(nx/2,10,max(nz/2,1)),tp, phase_field_new(nx/2,min(100,ny-1),max(nz/2,1))

4635  FORMAT(9(1X,e14.4))


end subroutine  equacao_campo_fase
