! ============================================================================
! karma3d.f90 -- modelo de campo de fase QUANTITATIVO para liga binaria diluida,
! crescimento isotermico livre (Karma, PRL 87, 115701, 2001; Echebarria, Folch,
! Karma & Plapp, PRE 70, 061604, 2004), em 3D, com GPU (OpenMP target offload).
!
! Variaveis adimensionais (comprimento W, tempo tau0):
!   phi in [-1,1]  (+1 solido, -1 liquido)
!   U = (c/c_l0 * 2/(1+k-(1-k)phi) - 1)/(1-k)   (supersaturacao adimensional)
!
!   a_s(n)^2 dphi/dt = div J + phi - phi^3 - lambda (1-phi^2)^2 U
!       J = a_s^2 grad(phi) + |grad phi|^2 a_s d(a_s)/d(grad phi)   (anisotropia cubica)
!       a_s(n) = 1 - 3 eps4 + 4 eps4 (nx^4+ny^4+nz^4)
!   A(phi) dU/dt = div[ D q(phi) grad U + a_t (1+(1-k)U) dphi/dt n ]
!                  + 1/2 (1+(1-k)U) dphi/dt
!       A = (1+k)/2 - (1-k)/2 phi,  q = (1-phi)/2 (sem difusao no solido),
!       a_t = 1/(2 sqrt 2) (corrente anti-trapping), n = grad(phi)/|grad(phi)|
!   Cinetica nula na interface: D = a2*lambda, d0 = a1*W/lambda,
!   a1 = 0.8839, a2 = 0.6267.
!
! Dimensional: W = (lambda/a1) d0, tau0 = a2 lambda W^2 / D_l,
!   d0 = Gamma/(|m|(1-k) c_l0), Gamma = sigma Tm / dH, c_l0 = (Tm-TP)/|m|,
!   Omega = (c_l0 - C0)/((1-k) c_l0), U(liquido longe) = -Omega.
!
! Dominio: 1/8 de semente no canto (i,j,k = 1), planos de simetria (espelho)
! em todas as faces. Celulas 1..N, fantasmas 0 e N+1.
! Saida: Tecplot POINT com X,Y,Z [m], PhaseField=(1+phi)/2, Temp=TP, Conc [fr. molar]
! -- as mesmas 6 variaveis do phasefield3d.
! ============================================================================
module karma_kernels
implicit none
contains

! Componente normal (direcao g1) do fluxo anisotropico J numa face, dado o
! gradiente g = (g1,g2,g3) na face.
pure function fluxn(g1, g2, g3, eps4) result(jn)
!$omp declare target
double precision, intent(in) :: g1, g2, g3, eps4
double precision :: jn, g2s, s4, as
g2s = g1*g1 + g2*g2 + g3*g3
if (g2s < 1.d-20) then
   jn = (1.d0-3.d0*eps4)**2 * g1
   return
end if
s4 = g1**4 + g2**4 + g3**4
as = 1.d0 - 3.d0*eps4 + 4.d0*eps4*s4/(g2s*g2s)
jn = as*( as*g1 + 16.d0*eps4*(g1**3 - g1*s4/g2s)/g2s )
end function fluxn

! Ruido uniforme em [-1,1) por (i,j,k,it), sem estado (seguro em GPU).
pure function ruido(i, j, k, it) result(r)
!$omp declare target
integer, intent(in) :: i, j, k, it
double precision :: r
integer(kind=8) :: h
integer(kind=8), parameter :: C1 = int(z'ff51afd7ed558ccd',8)
integer(kind=8), parameter :: C2 = int(z'c4ceb9fe1a85ec53',8)
integer(kind=8), parameter :: MASK63 = int(z'7fffffffffffffff',8)
h = int(i,8)*2654435761_8 + int(j,8)*2246822519_8 + int(k,8)*3266489917_8 &
  + int(it,8)*668265263_8 + 1_8
h = ieor(h, ishft(h,-33)); h = h*C1
h = ieor(h, ishft(h,-33)); h = h*C2
h = ieor(h, ishft(h,-33))
r = 2.d0*real(iand(h,MASK63),8)/real(MASK63,8) - 1.d0
end function ruido

end module karma_kernels


program karma3d
use karma_kernels
use omp_lib
implicit none

! ---- parametros (namelist karma.nml) ---------------------------------------
integer :: N = 200                 ! celulas por direcao
double precision :: lambda = 4.d0  ! acoplamento (W/d0 = lambda/a1)
double precision :: eps4 = 0.02d0  ! anisotropia cubica da energia de interface
double precision :: dx = 0.8d0     ! espacamento [W]
double precision :: dt = 0.008d0   ! passo de tempo [tau0]
double precision :: R0 = 8.d0      ! raio da semente [W]
double precision :: Fnoise = 0.d0  ! amplitude do ruido na eq. de fase
integer :: nsteps = 100000         ! passos maximos
integer :: nout = 2000             ! passos entre arquivos
double precision :: tip_stop = 0.9d0 ! para quando a ponta passar desta fracao do dominio
! material (Al-Cu, mesmos dados do entrada.dat)
double precision :: k = 0.14d0, Tm = 933.3d0, mliq = 640.d0, C0 = 0.0196d0
double precision :: TP = 909.5d0, sigma = 0.093d0, dH = 0.891002d9, Dl = 3.0d-9
namelist /karma/ N, lambda, eps4, dx, dt, R0, Fnoise, nsteps, nout, tip_stop, &
                 k, Tm, mliq, C0, TP, sigma, dH, Dl

double precision, parameter :: a1 = 0.8839d0, a2 = 0.6267d0
double precision, parameter :: at = 0.35355339059327373d0   ! 1/(2 sqrt 2)

double precision, allocatable :: phi(:,:,:), phin(:,:,:), U(:,:,:), Un(:,:,:), dphi(:,:,:)
double precision :: D, cl0, Omega, Gam, d0, W, tau0
double precision :: idx, idx2, x, y, z, r, g1, g2, g3, gm, as2, s4, divJ, pp
double precision :: Jxp, Jxm, Jyp, Jym, Jzp, Jzm, Fxp, Fxm, Fyp, Fym, Fzp, Fzm
double precision :: qf, cf, Au, rhsU, t0, t1
integer :: i, j, kk, it, ios, iu, tip_x, tip_d, nfile
logical :: use_gpu
character(len=64) :: fname

open(newunit=iu, file='karma.nml', status='old', iostat=ios)
if (ios == 0) then
   read(iu, nml=karma, iostat=ios)
   if (ios /= 0) stop 'erro lendo karma.nml'
   close(iu)
end if

D     = a2*lambda
cl0   = (Tm-TP)/mliq
Omega = (cl0-C0)/((1.d0-k)*cl0)
Gam   = sigma*Tm/dH
d0    = Gam/(mliq*(1.d0-k)*cl0)
W     = lambda/a1*d0
tau0  = a2*lambda*W*W/Dl

write(*,'(A)') '### karma3d -- campo de fase quantitativo (Karma 2001 / EFKP 2004) ###'
write(*,'(A,I0,A,F6.3,A,F8.3,A)') 'N = ', N, '^3,  dx = ', dx, ' W,  dominio = ', N*dx*W*1.d6, ' um'
write(*,'(A,F6.3,A,F6.3,A,F6.3,A,F7.3)') 'k = ', k, '  eps4 = ', eps4, '  lambda = ', lambda, '  D = ', D
write(*,'(A,F8.2,A,F8.5,A,F6.3)') 'TP = ', TP, ' K   c_l0 = ', cl0, '   Omega = ', Omega
write(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A)') 'd0 = ', d0, ' m   W = ', W, ' m   tau0 = ', tau0, ' s'
write(*,'(A,ES10.3,A,I0,A,F7.4,A)') 'dt = ', dt*tau0, ' s   (', nsteps, ' passos max)   Fnoise = ', Fnoise, ''
if (D*dt/(dx*dx) > 0.3d0/3.d0 .or. dt/(dx*dx) > 0.4d0/3.d0) &
   write(*,'(A)') 'AVISO: dt perto do limite de estabilidade explicita'

allocate(phi(0:N+1,0:N+1,0:N+1), phin(0:N+1,0:N+1,0:N+1), U(0:N+1,0:N+1,0:N+1), &
         Un(0:N+1,0:N+1,0:N+1), dphi(0:N+1,0:N+1,0:N+1))

! semente: 1/8 de esfera centrada no canto (plano de simetria em x=y=z=0)
!$omp parallel do collapse(2) private(i,j,kk,x,y,z,r)
do kk = 0, N+1
  do j = 0, N+1
    do i = 0, N+1
      x = (i-0.5d0)*dx ; y = (j-0.5d0)*dx ; z = (kk-0.5d0)*dx
      r = sqrt(x*x+y*y+z*z)
      phi(i,j,kk)  = -tanh((r-R0)/sqrt(2.d0))
      U(i,j,kk)    = -Omega
      phin(i,j,kk) = phi(i,j,kk)
      Un(i,j,kk)   = U(i,j,kk)
      dphi(i,j,kk) = 0.d0
    end do
  end do
end do

use_gpu = omp_get_num_devices() > 0
if (use_gpu) then
   write(*,'(A)') '### GPU offload ENABLED ###'
else
   write(*,'(A,I0,A)') '### CPU OpenMP (', omp_get_max_threads(), ' threads) ###'
end if

idx = 1.d0/dx ; idx2 = idx*idx
nfile = 0
t0 = omp_get_wtime()

!$omp target data map(tofrom: phi, U) map(alloc: phin, Un, dphi) if(use_gpu)

do it = 1, nsteps

  ! ---- fantasmas: espelho (fluxo nulo / simetria) em phi e U ----------------
  call espelho(phi) ; call espelho(U)

  ! ---- equacao de fase -------------------------------------------------------
  !$omp target teams distribute parallel do collapse(3) if(use_gpu) &
  !$omp private(i,j,kk,pp,g1,g2,g3,gm,s4,as2,Jxp,Jxm,Jyp,Jym,Jzp,Jzm,divJ)
  do kk = 1, N
    do j = 1, N
      do i = 1, N
        pp = phi(i,j,kk)
        ! fluxo nas 6 faces (gradiente na face: normal por diferenca, tangenciais pela media)
        Jxp = fluxn((phi(i+1,j,kk)-pp)*idx, &
                    0.25d0*idx*(phi(i,j+1,kk)+phi(i+1,j+1,kk)-phi(i,j-1,kk)-phi(i+1,j-1,kk)), &
                    0.25d0*idx*(phi(i,j,kk+1)+phi(i+1,j,kk+1)-phi(i,j,kk-1)-phi(i+1,j,kk-1)), eps4)
        Jxm = fluxn((pp-phi(i-1,j,kk))*idx, &
                    0.25d0*idx*(phi(i,j+1,kk)+phi(i-1,j+1,kk)-phi(i,j-1,kk)-phi(i-1,j-1,kk)), &
                    0.25d0*idx*(phi(i,j,kk+1)+phi(i-1,j,kk+1)-phi(i,j,kk-1)-phi(i-1,j,kk-1)), eps4)
        Jyp = fluxn((phi(i,j+1,kk)-pp)*idx, &
                    0.25d0*idx*(phi(i+1,j,kk)+phi(i+1,j+1,kk)-phi(i-1,j,kk)-phi(i-1,j+1,kk)), &
                    0.25d0*idx*(phi(i,j,kk+1)+phi(i,j+1,kk+1)-phi(i,j,kk-1)-phi(i,j+1,kk-1)), eps4)
        Jym = fluxn((pp-phi(i,j-1,kk))*idx, &
                    0.25d0*idx*(phi(i+1,j,kk)+phi(i+1,j-1,kk)-phi(i-1,j,kk)-phi(i-1,j-1,kk)), &
                    0.25d0*idx*(phi(i,j,kk+1)+phi(i,j-1,kk+1)-phi(i,j,kk-1)-phi(i,j-1,kk-1)), eps4)
        Jzp = fluxn((phi(i,j,kk+1)-pp)*idx, &
                    0.25d0*idx*(phi(i+1,j,kk)+phi(i+1,j,kk+1)-phi(i-1,j,kk)-phi(i-1,j,kk+1)), &
                    0.25d0*idx*(phi(i,j+1,kk)+phi(i,j+1,kk+1)-phi(i,j-1,kk)-phi(i,j-1,kk+1)), eps4)
        Jzm = fluxn((pp-phi(i,j,kk-1))*idx, &
                    0.25d0*idx*(phi(i+1,j,kk)+phi(i+1,j,kk-1)-phi(i-1,j,kk)-phi(i-1,j,kk-1)), &
                    0.25d0*idx*(phi(i,j+1,kk)+phi(i,j+1,kk-1)-phi(i,j-1,kk)-phi(i,j-1,kk-1)), eps4)
        divJ = (Jxp-Jxm+Jyp-Jym+Jzp-Jzm)*idx
        ! a_s(n)^2 no centro da celula (gradiente central)
        g1 = 0.5d0*idx*(phi(i+1,j,kk)-phi(i-1,j,kk))
        g2 = 0.5d0*idx*(phi(i,j+1,kk)-phi(i,j-1,kk))
        g3 = 0.5d0*idx*(phi(i,j,kk+1)-phi(i,j,kk-1))
        gm = g1*g1+g2*g2+g3*g3
        if (gm > 1.d-20) then
           s4  = (g1**4+g2**4+g3**4)/(gm*gm)
           as2 = (1.d0-3.d0*eps4+4.d0*eps4*s4)**2
        else
           as2 = (1.d0-3.d0*eps4)**2
        end if
        dphi(i,j,kk) = (divJ + pp - pp**3 - lambda*(1.d0-pp*pp)**2*U(i,j,kk))/as2 &
                     + Fnoise*ruido(i,j,kk,it)*(1.d0-pp*pp)/sqrt(dt)
        phin(i,j,kk) = pp + dt*dphi(i,j,kk)
      end do
    end do
  end do

  call espelho(dphi)

  ! ---- equacao de U (soluto) com anti-trapping ------------------------------
  !$omp target teams distribute parallel do collapse(3) if(use_gpu) &
  !$omp private(i,j,kk,pp,qf,cf,g1,g2,g3,gm,Fxp,Fxm,Fyp,Fym,Fzp,Fzm,Au,rhsU)
  do kk = 1, N
    do j = 1, N
      do i = 1, N
        pp = phi(i,j,kk)
        ! face x+
        g1 = (phi(i+1,j,kk)-pp)*idx
        g2 = 0.25d0*idx*(phi(i,j+1,kk)+phi(i+1,j+1,kk)-phi(i,j-1,kk)-phi(i+1,j-1,kk))
        g3 = 0.25d0*idx*(phi(i,j,kk+1)+phi(i+1,j,kk+1)-phi(i,j,kk-1)-phi(i+1,j,kk-1))
        gm = sqrt(g1*g1+g2*g2+g3*g3)
        qf = 0.25d0*(2.d0-pp-phi(i+1,j,kk))
        cf = 0.5d0*((1.d0+(1.d0-k)*U(i,j,kk))*dphi(i,j,kk)+(1.d0+(1.d0-k)*U(i+1,j,kk))*dphi(i+1,j,kk))
        Fxp = D*qf*(U(i+1,j,kk)-U(i,j,kk))*idx
        if (gm > 1.d-10) Fxp = Fxp + at*cf*g1/gm
        ! face x-
        g1 = (pp-phi(i-1,j,kk))*idx
        g2 = 0.25d0*idx*(phi(i,j+1,kk)+phi(i-1,j+1,kk)-phi(i,j-1,kk)-phi(i-1,j-1,kk))
        g3 = 0.25d0*idx*(phi(i,j,kk+1)+phi(i-1,j,kk+1)-phi(i,j,kk-1)-phi(i-1,j,kk-1))
        gm = sqrt(g1*g1+g2*g2+g3*g3)
        qf = 0.25d0*(2.d0-pp-phi(i-1,j,kk))
        cf = 0.5d0*((1.d0+(1.d0-k)*U(i,j,kk))*dphi(i,j,kk)+(1.d0+(1.d0-k)*U(i-1,j,kk))*dphi(i-1,j,kk))
        Fxm = D*qf*(U(i,j,kk)-U(i-1,j,kk))*idx
        if (gm > 1.d-10) Fxm = Fxm + at*cf*g1/gm
        ! face y+
        g2 = (phi(i,j+1,kk)-pp)*idx
        g1 = 0.25d0*idx*(phi(i+1,j,kk)+phi(i+1,j+1,kk)-phi(i-1,j,kk)-phi(i-1,j+1,kk))
        g3 = 0.25d0*idx*(phi(i,j,kk+1)+phi(i,j+1,kk+1)-phi(i,j,kk-1)-phi(i,j+1,kk-1))
        gm = sqrt(g1*g1+g2*g2+g3*g3)
        qf = 0.25d0*(2.d0-pp-phi(i,j+1,kk))
        cf = 0.5d0*((1.d0+(1.d0-k)*U(i,j,kk))*dphi(i,j,kk)+(1.d0+(1.d0-k)*U(i,j+1,kk))*dphi(i,j+1,kk))
        Fyp = D*qf*(U(i,j+1,kk)-U(i,j,kk))*idx
        if (gm > 1.d-10) Fyp = Fyp + at*cf*g2/gm
        ! face y-
        g2 = (pp-phi(i,j-1,kk))*idx
        g1 = 0.25d0*idx*(phi(i+1,j,kk)+phi(i+1,j-1,kk)-phi(i-1,j,kk)-phi(i-1,j-1,kk))
        g3 = 0.25d0*idx*(phi(i,j,kk+1)+phi(i,j-1,kk+1)-phi(i,j,kk-1)-phi(i,j-1,kk-1))
        gm = sqrt(g1*g1+g2*g2+g3*g3)
        qf = 0.25d0*(2.d0-pp-phi(i,j-1,kk))
        cf = 0.5d0*((1.d0+(1.d0-k)*U(i,j,kk))*dphi(i,j,kk)+(1.d0+(1.d0-k)*U(i,j-1,kk))*dphi(i,j-1,kk))
        Fym = D*qf*(U(i,j,kk)-U(i,j-1,kk))*idx
        if (gm > 1.d-10) Fym = Fym + at*cf*g2/gm
        ! face z+
        g3 = (phi(i,j,kk+1)-pp)*idx
        g1 = 0.25d0*idx*(phi(i+1,j,kk)+phi(i+1,j,kk+1)-phi(i-1,j,kk)-phi(i-1,j,kk+1))
        g2 = 0.25d0*idx*(phi(i,j+1,kk)+phi(i,j+1,kk+1)-phi(i,j-1,kk)-phi(i,j-1,kk+1))
        gm = sqrt(g1*g1+g2*g2+g3*g3)
        qf = 0.25d0*(2.d0-pp-phi(i,j,kk+1))
        cf = 0.5d0*((1.d0+(1.d0-k)*U(i,j,kk))*dphi(i,j,kk)+(1.d0+(1.d0-k)*U(i,j,kk+1))*dphi(i,j,kk+1))
        Fzp = D*qf*(U(i,j,kk+1)-U(i,j,kk))*idx
        if (gm > 1.d-10) Fzp = Fzp + at*cf*g3/gm
        ! face z-
        g3 = (pp-phi(i,j,kk-1))*idx
        g1 = 0.25d0*idx*(phi(i+1,j,kk)+phi(i+1,j,kk-1)-phi(i-1,j,kk)-phi(i-1,j,kk-1))
        g2 = 0.25d0*idx*(phi(i,j+1,kk)+phi(i,j+1,kk-1)-phi(i,j-1,kk)-phi(i,j-1,kk-1))
        gm = sqrt(g1*g1+g2*g2+g3*g3)
        qf = 0.25d0*(2.d0-pp-phi(i,j,kk-1))
        cf = 0.5d0*((1.d0+(1.d0-k)*U(i,j,kk))*dphi(i,j,kk)+(1.d0+(1.d0-k)*U(i,j,kk-1))*dphi(i,j,kk-1))
        Fzm = D*qf*(U(i,j,kk)-U(i,j,kk-1))*idx
        if (gm > 1.d-10) Fzm = Fzm + at*cf*g3/gm

        Au   = 0.5d0*(1.d0+k) - 0.5d0*(1.d0-k)*pp
        rhsU = (Fxp-Fxm+Fyp-Fym+Fzp-Fzm)*idx + 0.5d0*(1.d0+(1.d0-k)*U(i,j,kk))*dphi(i,j,kk)
        Un(i,j,kk) = U(i,j,kk) + dt*rhsU/Au
      end do
    end do
  end do

  ! ---- avanca no tempo -------------------------------------------------------
  !$omp target teams distribute parallel do collapse(3) if(use_gpu) private(i,j,kk)
  do kk = 1, N
    do j = 1, N
      do i = 1, N
        phi(i,j,kk) = max(-1.d0, min(1.d0, phin(i,j,kk)))
        U(i,j,kk)   = Un(i,j,kk)
      end do
    end do
  end do

  ! ---- ponta: maior i com phi>0 no eixo x (j=k=1) e na diagonal -------------
  if (mod(it, nout) == 0 .or. it == nsteps) then
     !$omp target update from(phi, U) if(use_gpu)
     tip_x = 0 ; tip_d = 0
     do i = 1, N
        if (phi(i,1,1) > 0.d0) tip_x = i
        if (phi(i,i,i) > 0.d0) tip_d = i
     end do
     t1 = omp_get_wtime()
     write(*,'(A,I8,A,F10.1,A,ES10.3,A,F7.1,A,F7.1,A,F6.3,A,F8.1,A)') 'it=', it, '  t=', it*dt, &
          ' tau0 (', it*dt*tau0, ' s)  ponta_eixo=', tip_x*dx, ' W  ponta_diag=', &
          tip_d*dx*sqrt(3.d0), ' W  razao=', tip_x/max(tip_d*sqrt(3.d0),1.d0), '  wall=', t1-t0, ' s'
     call flush(6)
     nfile = nfile + 1
     write(fname,'(A,I7.7,A)') 'output/karma_', it, '.dat'
     call escreve(fname)
     if (tip_x > tip_stop*N) then
        write(*,'(A)') '### ponta chegou perto da parede -- fim ###'
        exit
     end if
  end if
end do

!$omp end target data

contains

subroutine espelho(f)
double precision, intent(inout) :: f(0:N+1,0:N+1,0:N+1)
integer :: a, b
!$omp target teams distribute parallel do collapse(2) if(use_gpu) private(a,b)
do b = 0, N+1
  do a = 0, N+1
    f(0,a,b)   = f(1,a,b)   ; f(N+1,a,b) = f(N,a,b)
  end do
end do
!$omp target teams distribute parallel do collapse(2) if(use_gpu) private(a,b)
do b = 0, N+1
  do a = 0, N+1
    f(a,0,b)   = f(a,1,b)   ; f(a,N+1,b) = f(a,N,b)
  end do
end do
!$omp target teams distribute parallel do collapse(2) if(use_gpu) private(a,b)
do b = 0, N+1
  do a = 0, N+1
    f(a,b,0)   = f(a,b,1)   ; f(a,b,N+1) = f(a,b,N)
  end do
end do
end subroutine espelho

subroutine escreve(nome)
character(len=*), intent(in) :: nome
integer :: iou, a, b, c
double precision :: cc
open(newunit=iou, file=trim(nome), status='replace', action='write')
write(iou,'(A)') 'TITLE="Phase Field (Karma)"'
write(iou,'(A)') 'variables= "X","Y","Z","PhaseField","Temp","Conc_Solute"'
write(iou,'(A)') 'ZONE'
write(iou,'(A,I0,A,I0,A,I0,A)') ' T=" " ,I= ', N, ',J= ', N, ',K= ', N, ',F=POINT'
do c = 1, N
  do b = 1, N
    do a = 1, N
      cc = cl0*0.5d0*(1.d0+k-(1.d0-k)*phi(a,b,c))*(1.d0+(1.d0-k)*U(a,b,c))
      write(iou,'(6(1X,ES15.7E3))') (a-1)*dx*W, (b-1)*dx*W, (c-1)*dx*W, &
                                    0.5d0*(1.d0+phi(a,b,c)), TP, cc
    end do
  end do
end do
close(iou)
end subroutine escreve

end program karma3d
