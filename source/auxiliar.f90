SUBROUTINE AUXILIAR
Use Variaveis
IMPLICIT NONE
Double precision::Phi_i, h_phi, numerador, denominador,function_1, Function_Integer, Function_Integer_Old, Integral_Trapezio_Repetida, Qsi, Inv_M0
Double precision:: E00, W00, f_a, f_b

!===============dados de entrada=============================================



!========================CALCULO DO TRAPÉZIO REPETIDA=====================================================
Function_Integer_Old = 0.d0

  phi_i = 0.10d0

10  h_phi = phi_i**3*(10.d0-15.d0*phi_i+6.d0*phi_i*phi_i)   ! quintico (anterior: phi_i*phi_i*(3.d0-2.d0*phi_i))
    function_1 = phi_i*(1.d0-phi_i)

numerador = h_phi*(1.d0-h_phi)
denominador = (  (1.d0-h_phi)*Cl_eq*(1.d0-Cl_eq) + h_Phi*Cs_eq*(1.d0-Cs_eq)  )*function_1

Function_Integer = (numerador/denominador)
Function_Integer_Old = Function_Integer_Old + Function_Integer
if (abs(phi_i-0.10d0)<1.d-6) f_a = Function_Integer
f_b = Function_Integer

phi_i = phi_i+0.01d0

IF (Phi_i<=0.9001d0) then
go to 10
end if

! Trapezio de 0.10 a 0.90 com os extremos calculados (antes: 54.39180136 e 234.073053 fixos,
! validos so para um TP e o h cubico)
Integral_Trapezio_Repetida = 0.01d0*( Function_Integer_Old - 0.5d0*(f_a+f_b) )

!=======================FIM DO CALCULO DA INTEGRAL========================================================


Qsi = R*TP*(1.d0/Vm)*(Cl_eq-Cs_eq)*(Cl_eq-Cs_eq)*Integral_Trapezio_Repetida


Inv_M0 = (E0*E0/Sigma) * ( (Beta*R*Tp*(1-Ke)/(Vm*me)) + (E0*Qsi/(D_liq*sqrt(2.d0*W0)))  )


M0 = ajuste*1.d0/Inv_M0


! M0 = (D_liq/(E0*E0))*0.1D-00



END SUBROUTINE AUXILIAR
