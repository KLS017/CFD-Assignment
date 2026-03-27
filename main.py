import numpy as np
import matplotlib.pyplot as plt

Lx = Ly = 1.0
Nx = Ny = 50
dx = Lx / Nx
dy = Ly / Ny
nframes = 40000
dt = 0.001
ipause = 0.1

x = np.linspace(0.5 * dx, Lx - 0.5 * dx, Nx)
y = np.linspace(0.5 * dy, Ly - 0.5 * dy, Ny)
X, Y = np.meshgrid(x, y) # Cell-center coordinates associated with I(x, y)
#inlet and outlet locations of Q1.4
ymin_in, ymax_in = 0.86, 0.94
ymin_out, ymax_out = 0.06, 0.14
#alpha values
alpha_liquid = 1e-1
alpha_wall = 1e-8

I = np.load("maze_geometry.npy")

c = np.zeros((Ny+2, Nx+2))
alpha = np.zeros((Ny+2, Nx+2))

#determine alphas for the walls
alpha[1:-1, 1:-1] = alpha_liquid*(1-I) + alpha_wall * I
alpha[0, :] = alpha[1, :]
alpha[-1, :] = alpha[-2, :]
alpha[:, 0] = alpha[:, 1]
alpha[:, -1] = alpha[:, -2]


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



#frame rate = nframes * ipause


for n in range(nframes):
    for j in range(1, Ny + 1):
        yj = (j - 0.5) * dy
        if yj >= ymin_in and yj <= ymax_in:
            c[j, 0] = 2 - c[j, 1]
        else:
            c[j, 0] = c[j,1]
        if yj >= ymin_out and yj <= ymax_out:
            c[j, Nx + 1] = -c[j, Nx]
        else:
            c[j, Nx + 1] = c[j, Nx]
    c[0, :] = c[1, :]
    c[Ny + 1, :] = c[Ny, :]

    
    c_new = c.copy()
    for j in range(1, Ny + 1):
        for i in range(1, Nx + 1):
            alpha_up    = 2 * alpha[j, i] * alpha[j+1, i] / (alpha[j, i] + alpha[j+1, i])
            alpha_down  = 2 * alpha[j, i] * alpha[j-1, i] / (alpha[j, i] + alpha[j-1, i])
            alpha_right = 2 * alpha[j, i] * alpha[j, i+1] / (alpha[j, i] + alpha[j, i+1])
            alpha_left  = 2 * alpha[j, i] * alpha[j, i-1] / (alpha[j, i] + alpha[j, i-1])
            
            c_new[j,i] = c[j,i] + dt * (
                (alpha_right * (c[j,i+1] - c[j,i]) - alpha_left * (c[j,i] - c[j,i-1])) / dx**2
                +
                (alpha_up * (c[j+1,i] - c[j,i]) - alpha_down * (c[j,i] - c[j-1,i])) / dy**2
            )


    c = c_new



        # inlet & outlet coefficients of diffusion
        # c[1, 44] = c[1, 45] = c[1, 46] = c[1, 47] = c[1, 48] = 1
        # c[4,51] = c[5,51] =c[6,51] =c[7,51] =c[8,51] = 0

    
    if n == 1/dt or n == 5/dt or n == 15/dt:  
        plt.clf()
        maze_walls = np.ma.masked_where(I == 0, I)
        plt.imshow(
            #no ghost cell plotting
            c[1:-1, 1:-1],
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