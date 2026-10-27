module ruido_mod
! Ruido espacialmente suave para a equacao de fase: valores sorteados (hash
! inteiro, sem estado) numa grade grossa de passo LR celulas, renovados a cada
! TR passos, e interpolados bilinearmente. Perturba na escala dos bracos
! secundarios (~um) sem enrugar a interface na escala da malha.
implicit none
integer, parameter :: LR = 6, TR = 20
contains

pure function hash01(a, b, c) result(r)
integer, intent(in) :: a, b, c
double precision :: r
integer(kind=8) :: h
integer(kind=8), parameter :: C1 = int(z'ff51afd7ed558ccd',8)
integer(kind=8), parameter :: C2 = int(z'c4ceb9fe1a85ec53',8)
integer(kind=8), parameter :: MASK63 = int(z'7fffffffffffffff',8)
h = int(a,8)*2654435761_8 + int(b,8)*2246822519_8 + int(c,8)*3266489917_8 + 1_8
h = ieor(h, ishft(h,-33)); h = h*C1
h = ieor(h, ishft(h,-33)); h = h*C2
h = ieor(h, ishft(h,-33))
r = real(iand(h,MASK63),8)/real(MASK63,8)
end function hash01

! Valor em [0,1) na celula (i,j) no passo it.
pure function ruido_suave(i, j, it) result(r)
integer, intent(in) :: i, j, it
double precision :: r, fx, fy
integer :: i0, j0, tb
i0 = i/LR ; j0 = j/LR ; tb = it/TR
fx = dble(mod(i,LR))/dble(LR) ; fy = dble(mod(j,LR))/dble(LR)
r = (1.d0-fx)*(1.d0-fy)*hash01(i0,  j0,  tb) + fx*(1.d0-fy)*hash01(i0+1,j0,  tb) &
  + (1.d0-fx)*fy       *hash01(i0,  j0+1,tb) + fx*fy       *hash01(i0+1,j0+1,tb)
end function ruido_suave

end module ruido_mod
