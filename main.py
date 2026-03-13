import numpy as np
import matplotlib.pyplot as plt

Lx = Ly = 1.0
Nx = Ny = 50
dx = Lx / Nx
dy = Ly / Ny
x = np.linspace(0.5 * dx, Lx - 0.5 * dx, Nx)
y = np.linspace(0.5 * dy, Ly - 0.5 * dy, Ny)
X, Y = np.meshgrid(x, y) # Cell-center coordinates associated with I(x, y)

I = np.load("maze_geometry.npy")

im = plt.imshow(
    I,
    origin="lower",
    extent=[0.0, Lx, 0.0, Ly],
    cmap="binary",
    vmin=0,
    vmax=1,
    interpolation="nearest",
    aspect="equal",
)

plt.colorbar(im, ticks=[0, 1], label="$I(x,y)$")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.show()


