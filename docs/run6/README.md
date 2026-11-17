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

Passo a passo (WSL, GPU NVIDIA):

```bash
cp -r <repo> ~/phasefield && cd ~/phasefield && bash phasefield_gpu.sh build
mkdir -p ~/run6/output && cd ~/run6
cp ~/phasefield/phasefield3d_gpu_wsl ~/phasefield/phasefield_gpu.sh .
cp ~/phasefield/docs/run6/entrada_run6.dat entrada.dat
bash phasefield_gpu.sh run      # 0 -> 52,47 ms (~2 dias na RTX 4050)
# run6b: nova pasta com entrada_run6b.dat e restart.dat -> saida de 52,5 ms do run6
# run6c: nova pasta com entrada_run6c.dat e restart.dat -> saida 22501 do run6b
```

## Figuras

Geradas com as macros de `tecplot/` (detalhes no README principal, secao
"Figuras 3D no Tecplot"):

- `tecplot/dendrita_wbm_fase.mcr`: solido vermelho (iso-superficie phi = 0,5, cor solida
  com iluminacao Gouraud) e paredes x = 0, y = 0 coloridas pelo PhaseField (azul = liquido).
- `tecplot/dendrita_wbm.mcr`: a mesma iso-superficie colorida pela concentracao de Cu.

Copie o arquivo de saida para uma pasta do Windows, ajuste `|ARQ|` e `|PNG|` no topo da
macro e rode `tec360.exe -b -p <macro>`.

## Resultados

| t (ms) | Figura | Altura da ponta (phi >= 0,5) |
|---|---|---|
| 30 | `fase_30ms.png` | 43 um: colunas em ogiva |
| 52,5 | `fase_52.5ms.png` | 99 um: bracos secundarios nitidos |
| 75 | `fase_75ms.png`, `Cu_75ms.png` | 158 um (90,9 + 67,5): bracos em degraus |

Velocidade da ponta ~2,5 um/ms, cerca de 1,8x a do 2D (`caso2D_75ms.png`: ponta a ~87 um
aos 75 ms). Superficie lisa, sem enrugamento. Os furos nas colunas cortadas pela parede
sao bolsoes de liquido interdendritico.
