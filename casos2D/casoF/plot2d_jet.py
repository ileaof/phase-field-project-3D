"""Cu de um arquivo 2D do phasefield em colormap jet (estilo das figuras do usuario)."""
import re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
num = lambda s: float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", s))
with open(sys.argv[1]) as f:
    for _ in range(3): f.readline()
    z = f.readline(); nx = int(z.split("I=")[1].split(",")[0]); ny = int(z.split("J=")[1].split(",")[0])
    d = np.array([[num(t) for t in f.readline().split()] for _ in range(nx*ny)])
cu = d[:, 4].reshape(ny, nx) * 100
x, y = d[:nx, 0]*1e6, d[::nx, 1]*1e6
fig, ax = plt.subplots(figsize=(6, 7), constrained_layout=True)
im = ax.imshow(cu, origin="lower", extent=[x[0], x[-1], y[0], y[-1]], cmap="jet")
fig.colorbar(im, ax=ax, shrink=0.7, label="Cu (% molar)")
ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)"); ax.set_title(sys.argv[3] if len(sys.argv) > 3 else "")
fig.savefig(sys.argv[2], dpi=160); print("salvo", sys.argv[2])
