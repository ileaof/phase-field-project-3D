# Caso F (2D) — dendritas colunares com o código 2D do usuário

Código 2D do projeto `Visual Studio 2022\Projects\PhaseFieldAlCu` com três
alterações (o resto é idêntico ao original):

| Arquivo | Alteração |
|---|---|
| `equacao_campo_fase.f90` | h(φ) quíntico; ruído de média zero `(2*ran-1)`; bloco `TX` desligado (estourava o array) |
| `concentracao_soluto.f90` | h(φ) quíntico |
| `temperatura.f90` | h'(φ) quíntico no termo fonte |

Parâmetros (`entrada.dat`; valor usado primeiro, anteriores ao lado): malha 300×400,
dx = dy = 1,5e-7 m (45×60 µm); T0 = 905 K (o resfriamento de 19,84 K/s do código leva
a 904,2 K em 40 ms); AAA = 0,01; DELTA_EPSLON = 0,03; ajuste = 2; Thin_Interface = 4;
dt = dx²/(5·D_liq) = 1,5e-6 s; stoptime = 0,04 s; saída a cada 2000 passos;
frente plana em j ≤ 4 crescendo em +y.

Resultado: `casoF_40ms_Cu.png` — dendritas colunares de ponta aguda, competição
entre troncos e Cu segregado no líquido interdendrítico.

## Rodar (Windows, Intel ifx)

```bat
call ..\..\setvars.bat
mkdir output
ifx -O2 source\variaveis.f90 source\W0_E0.f90 source\auxiliar.f90 source\concentracao_soluto.f90 ^
    source\entrada.f90 source\equacao_campo_fase.f90 source\principal.f90 source\saida.f90 ^
    source\temperatura.f90 source\Calc_TemperatureNew.f90 -o pf2d.exe
pf2d.exe
python plot2d_jet.py "output\<arquivo>.dat" figura.png "titulo"
```
