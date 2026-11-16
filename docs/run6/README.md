# run6: dendritas colunares 3D com os valores default do 2D

Modelo WBM 3D (`source/`) com os mesmos parametros do codigo 2D original do usuario
(`casos2D/`): Al-Cu, T0 = 922 K, resfriamento continuo de 19,84 K/s, AAA = 0,03,
h(phi) cubico e ruido positivo celula a celula (linha 43 do entrada.dat: `0 0`).
Malha 340 x 340 x 700 (51 x 51 x 105 um), dx = 1,5e-7 m, frente plana na base crescendo em +z.
GPU RTX 4050 (6 GB, ~4,4 GB usados), nvfortran no WSL, ~3 s/passo, saida a cada 7500 passos (~7,5 ms).

## Como repetir

O ruido e deterministico (`pseudo_random(i,j,k,it)` em `variaveis.f90`): os mesmos fontes e
entrada.dat reproduzem a mesma simulacao. O dominio de 105 um so comporta as pontas ate ~52,5 ms;
depois a simulacao continua com o **referencial movel** (`tipo_inicial = 4`, segundo valor da
linha 41 = numero de camadas descartadas na base; o topo e preenchido com liquido).

| Etapa | entrada.dat | Condicao inicial | Tempo real | Deslocamento z |
|---|---|---|---|---|
| run6  | `entrada_run6.dat`  | frente plana (3) | 0 - 52,47 ms | 0 |
| run6b | `entrada_run6b.dat` | restart.dat = saida de 52,5 ms do run6, 450 camadas | 52,47 + t | +67,5 um |
| run6c | `entrada_run6c.dat` | restart.dat = saida de 22501 passos do run6b, 300 camadas | 74,975 + t | +112,5 um |

Em cada continuacao: T0 = 922 - 19,84 * (tempo real em s) e stoptime = 85 ms - tempo real.
Na retomada o contador de passos volta a 1 (e com ele a sequencia do ruido); para resultados
identicos, retome nos mesmos pontos.

## Resultados

| t (ms) | Figura | Altura da ponta (phi >= 0,5) |
|---|---|---|
| 30 | `fase_30ms.png` | 43 um: colunas em ogiva |
| 52,5 | `fase_52.5ms.png` | 99 um: bracos secundarios nitidos |
| 75 | `fase_75ms.png`, `Cu_75ms.png` | 158 um (90,9 + 67,5): bracos em degraus |

Velocidade da ponta ~2,5 um/ms, cerca de 1,8x a do 2D (`caso2D_75ms.png`: ponta a ~87 um
aos 75 ms). Superficie lisa, sem enrugamento. Os furos nas colunas cortadas pela parede
sao bolsoes de liquido interdendritico.
