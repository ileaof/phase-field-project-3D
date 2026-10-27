# Caso O (2D) — dendritas colunares com ruído suave

Código 2D do usuário (projeto `Visual Studio 2022\Projects\PhaseFieldAlCu`) com as
alterações do caso F (h(φ) quíntico, ruído de média zero, bloco `TX` desligado) e
**ruído espacialmente suave** no lugar de `call random(ran)`:

| Arquivo | Alteração em relação ao caso F |
|---|---|
| `ruido_suave.f90` (novo) | valores sorteados a cada 6 células, renovados a cada 20 passos e interpolados bilinearmente |
| `equacao_campo_fase.f90` | `ran = ruido_suave(i,j,it)` em vez de `call random(ran)`; `use ruido_mod` |

O ruído célula a célula enruga a interface na escala da malha (amplitude alta) ou não dispara
braços secundários (amplitude baixa). O ruído suave perturba na escala dos braços (~1 µm):
superfície lisa **e** início de braços secundários. Com o ruído suave a amplitude precisa ser
~10× menor que a do ruído célula a célula (AAA = 0,01 e 0,03 com ruído suave formaram "espuma").

Parâmetros (`entrada.dat`; valor usado primeiro, anteriores ao lado): **1200×800**,
dx = dy = 1,5e-7 m (180×120 µm); T0 = 905 K; **AAA = 0,001**; DELTA_EPSLON = 0,03; ajuste = 2;
dt = dx²/(5·D_liq) = 1,5e-6 s; stoptime = 0,08 s; saída a cada 4000 passos; frente plana em
j ≤ 4 crescendo em +y. O primeiro teste (300×800, 45 µm de largura) está em
`docs/run4_ruido_suave/caso2D_O_Cu.png`.

## Rodar (Windows, Intel ifx)

```bat
call ..\..\setvars.bat
mkdir output
ifx -O2 source\variaveis.f90 source\ruido_suave.f90 source\W0_E0.f90 source\auxiliar.f90 ^
    source\concentracao_soluto.f90 source\entrada.f90 source\equacao_campo_fase.f90 ^
    source\principal.f90 source\saida.f90 source\temperatura.f90 source\Calc_TemperatureNew.f90 -o pf2d.exe
pf2d.exe
python plot2d_jet.py "output\<arquivo>.dat" figura.png "titulo"
```
