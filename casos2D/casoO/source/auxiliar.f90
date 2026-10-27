SUBROUTINE AUXILIAR
Use Variaveis
IMPLICIT NONE
Double precision::Phi_i, h_phi, numerador, denominador,function_1, Function_Integer, Function_Integer_Old, Integral_Trapezio_Repetida, Qsi, Inv_M0
Double precision:: E00, W00

!===============dados de entrada=============================================



!========================CALCULO DO TRAPÉZIO REPETIDA=====================================================
Function_Integer_Old = 0.d0

  phi_i = 0.11d0

10  h_phi = phi_i*phi_i*(3.d0-2.d0*phi_i)
    function_1 = phi_i*(1.d0-phi_i)

numerador = h_phi*(1.d0-h_phi)
denominador = (  (1.d0-h_phi)*Cl_eq*(1.d0-Cl_eq) + h_Phi*Cs_eq*(1.d0-Cs_eq)  )*function_1

Function_Integer = (numerador/denominador)
Function_Integer_Old = Function_Integer_Old + Function_Integer

phi_i = phi_i+0.01d0

IF (Phi_i<=0.8999d0) then
go to 10
end if

Integral_Trapezio_Repetida = (0.01/2.d0)*( 54.39180136 + 2.d0*(Function_Integer_Old) + 234.073053d0 )

!=======================FIM DO CALCULO DA INTEGRAL========================================================


Qsi = R*TP*(1.d0/Vm)*(Cl_eq-Cs_eq)*(Cl_eq-Cs_eq)*Integral_Trapezio_Repetida


Inv_M0 = (E0*E0/Sigma) * ( (Beta*R*Tp*(1-Ke)/(Vm*me)) + (E0*Qsi/(D_liq*sqrt(2.d0*W0)))  )


M0 = ajuste*1.d0/Inv_M0


! M0 = (D_liq/(E0*E0))*0.1D-00



END SUBROUTINE AUXILIAR
