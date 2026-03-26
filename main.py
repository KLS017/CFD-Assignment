import numpy as np
import matplotlib.pyplot as plt

Lx = Ly = 1.0
Nx = Ny = 50
dx = Lx / Nx
dy = Ly / Ny
x = np.linspace(0.5 * dx, Lx - 0.5 * dx, Nx)
y = np.linspace(0.5 * dy, Ly - 0.5 * dy, Ny)
X, Y = np.meshgrid(x, y) # Cell-center coordinates associated with I(x, y)
#inlet and outlet locations of Q1.4
ymin_in, ymax_in = 0.875, 0.955
ymin_out, ymax_out = 0.06, 0.14

I = np.load("maze_geometry.npy")


# im = plt.imshow(
#     I,
#     origin="lower",
#     extent=[0.0, Lx, 0.0, Ly],
#     cmap="binary",
#     vmin=0,
#     vmax=1,
#     interpolation="nearest",
#     aspect="equal",
# )

# plt.colorbar(im, ticks=[0, 1], label="$I(x,y)$")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.tight_layout()
# plt.show()

c = np.zeros((Ny+2, Nx+2))



plt.figure()
plt.ion()

nframes = 200
ipause = 0.01
#frame rate = nframes * ipause


for n in range(nframes):
    
    for j in range(1, Ny + 1):
        
        yj = (j - 0.5) * dy
        if yj >= ymin_in and yj <= ymax_in:
            c[j, 0] = 1 - c[j, 1]
        else:
            c[j, 0] = c[j,1]
        if yj >= ymin_out and yj <= ymax_out:
            c[j, Nx + 1] = 0
        else:
            c[j, Nx + 1] = 0
        c[0, :] = c[1, :]
        c[Ny + 1, :] = c[Ny, :]
        

        

        # inlet & outlet coefficients of diffusion
        # c[1, 44] = c[1, 45] = c[1, 46] = c[1, 47] = c[1, 48] = 1
        # c[4,51] = c[5,51] =c[6,51] =c[7,51] =c[8,51] = 0

    
    if n % 50 == 0:  
        print(c[1, 44])
        plt.clf()
        maze_walls = np.ma.masked_where(I == 0, I)
        plt.imshow(
            c,
            origin="lower",
            extent=[0.0, Lx, 0.0, Ly],
            cmap="viridis",
            aspect="equal",
        )
        plt.imshow(
            maze_walls,
            origin="lower",
            extent=[0.0, Lx, 0.0, Ly],
            cmap="gray_r",
            alpha=0.15,
            interpolation="nearest",
            aspect="equal",
        )
        plt.colorbar()
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Transient concentration")
        plt.pause(ipause)
        plt.ioff()    
        plt.show() 