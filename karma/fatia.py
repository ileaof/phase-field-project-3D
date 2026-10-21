"""Imagem de um corte 2D (plano z = cte) de um arquivo do karma3d.

Le so as linhas do plano pedido (o arquivo Tecplot POINT tem i mais rapido,
depois j, depois k), espelha nos dois planos de simetria para mostrar a secao
inteira da dendrita e salva um PNG com a fase e a concentracao lado a lado.

Uso (no WSL, na pasta do run):
    python3 fatia.py output/karma_0045000.dat            # plano z = 0
    python3 fatia.py output/karma_0045000.dat --k 20     # plano k = 20
"""
import argparse
import itertools

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("arquivo")
ap.add_argument("--k", type=int, default=1, help="indice do plano z (1 = plano de simetria)")
ap.add_argument("--saida", default=None, help="PNG de saida (padrao: <arquivo>_k<k>.png)")
a = ap.parse_args()

with open(a.arquivo) as f:
    for _ in range(3):
        f.readline()
    zona = f.readline()                      # ' T=" " ,I= N,J= N,K= N,F=POINT'
    N = int(zona.split("I=")[1].split(",")[0])
    linhas = itertools.islice(f, (a.k - 1) * N * N, a.k * N * N)
    d = np.loadtxt(linhas).reshape(N, N, 6)  # [j, i, var]

x = d[0, :, 0] * 1e6                         # um
fase, conc = d[:, :, 3], d[:, :, 5]

def espelha(c):
    """Quadrante -> secao inteira (planos de simetria x=0 e y=0)."""
    c = np.concatenate([c[::-1, :], c], axis=0)
    return np.concatenate([c[:, ::-1], c], axis=1)

ext = [-x[-1], x[-1], -x[-1], x[-1]]
fig, ax = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
im0 = ax[0].imshow(espelha(fase), origin="lower", extent=ext, cmap="gray_r", vmin=0, vmax=1)
ax[0].set_title("Fase (1 = sólido)")
im1 = ax[1].imshow(espelha(conc) * 100, origin="lower", extent=ext, cmap="viridis")
ax[1].contour(espelha(fase), levels=[0.5], colors="w", linewidths=0.6, extent=ext)
ax[1].set_title("Cu (% molar), contorno φ = 0,5")
for axi in ax:
    axi.set_xlabel("x (µm)"); axi.set_ylabel("y (µm)")
fig.colorbar(im0, ax=ax[0], shrink=0.8)
fig.colorbar(im1, ax=ax[1], shrink=0.8)
saida = a.saida or a.arquivo.rsplit(".", 1)[0] + f"_k{a.k}.png"
fig.savefig(saida, dpi=200)
print("salvo:", saida)
