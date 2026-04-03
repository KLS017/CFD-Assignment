import numpy as np
import matplotlib.pyplot as plt

import time

Lx = Ly = 1.0
Nx = Ny = 50
dx = Lx / Nx
dy = Ly / Ny
t = 50
dt = 0.01 #should be a dt calculation here with use of the alphas but wasn't sure how
nframes = int(t / dt)
ipause = 0.1
tol = 1e-6
max_iter = 200
skip = 40
# X, Y = np.meshgrid(x, y) # Cell-center coordinates associated with I(x, y)
#inlet and outlet locations of Q1.4
ymin_in, ymax_in = 0.86, 0.94
ymin_out, ymax_out = 0.06, 0.14
#alpha values
alpha_liquid = 1e-2
alpha_wall = 1e-8


x = np.linspace(0.5 * dx, Lx - 0.5 * dx, Nx)
y = np.linspace(0.5 * dy, Ly - 0.5 * dy, Ny)
diagnostic_quantity = np.zeros((Ny, Nx))

I = np.load("maze_geometry.npy")

c = np.zeros((Ny+2, Nx+2))
alpha = np.zeros((Ny+2, Nx+2))

#determine alphas for the walls
alpha[1:-1, 1:-1] = alpha_liquid*(1-I) + alpha_wall * I
#determine alpha for the ghost cells
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

    

plt.ion()
fig, ax = plt.subplots()

start = time.time()
for n in range(nframes+1):
    c_new = c.copy()
    c_new[0, :] = c[1, :]
    c_new[Ny + 1, :] = c[Ny, :]  
    
    for it in range(max_iter):
        c_new[0, :] = c_new[1, :] 
        c_new[Ny + 1, :] = c_new[Ny, :]
        c_old_iter = c_new.copy()
        for j in range(1, Ny + 1):
            yj = (j - 0.5) * dy
            if yj >= ymin_in and yj <= ymax_in:
                c_new[j, 0] = 2 - c_new[j, 1]
            else:
                c_new[j, 0] = c_new[j,1]
            if yj >= ymin_out and yj <= ymax_out:
                c_new[j, Nx + 1] = -c_new[j, Nx]
            else:
                c_new[j, Nx + 1] = c_new[j, Nx]
            for i in range(1, Nx + 1):
                alpha_up    = 2 * alpha[j, i] * alpha[j+1, i] / (alpha[j, i] + alpha[j+1, i])
                alpha_down  = 2 * alpha[j, i] * alpha[j-1, i] / (alpha[j, i] + alpha[j-1, i])
                alpha_right = 2 * alpha[j, i] * alpha[j, i+1] / (alpha[j, i] + alpha[j, i+1])
                alpha_left  = 2 * alpha[j, i] * alpha[j, i-1] / (alpha[j, i] + alpha[j, i-1])
                Ar = dt * alpha_right / dx**2
                Al = dt * alpha_left  / dx**2
                Au = dt * alpha_up    / dy**2
                Ad = dt * alpha_down  / dy**2

                den = 1 + Ar + Al + Au + Ad
                c_new[j, i] = (
                    c[j, i]
                    + Al * c_new[j, i-1]
                    + Ad * c_new[j-1, i]
                    + Ar * c_old_iter[j, i+1]
                    + Au * c_old_iter[j+1, i]
                ) / den

        iter_err = np.linalg.norm(c_new[1:Ny+1, 1:Nx+1]- c_old_iter[1:Ny+1, 1:Nx+1], ord=2)  #DELETE LATER: added [1:Ny+1, 1:Nx+1] to ensure ghost cells aren't counted, added ord=2 cause pedro does that as well
        if iter_err < tol:
            break
    if abs(tn1 % 0.5) < 1e-12:
        print(f"t={tn1:5.2f} | GS={it:3d} | Δ={np.max(np.abs(c_new - c)):.3e}")

            # alpha_up    = 2 * alpha[j, i] * alpha[j+1, i] / (alpha[j, i] + alpha[j+1, i])
            # alpha_down  = 2 * alpha[j, i] * alpha[j-1, i] / (alpha[j, i] + alpha[j-1, i])
            # alpha_right = 2 * alpha[j, i] * alpha[j, i+1] / (alpha[j, i] + alpha[j, i+1])
            # alpha_left  = 2 * alpha[j, i] * alpha[j, i-1] / (alpha[j, i] + alpha[j, i-1])
            # c[0, :] = c[1, :]
            # c[Ny + 1, :] = c[Ny, :]  
            # c_new[j,i] = c[j,i] + dt * (
            #     (alpha_right * (c[j,i+1] - c[j,i]) - alpha_left * (c[j,i] - c[j,i-1])) / dx**2
            #     +
            #     (alpha_up * (c[j+1,i] - c[j,i]) - alpha_down * (c[j,i] - c[j-1,i])) / dy**2
            # )


    c = c_new
    
    # if n % skip == 0:
    #     c_copy = c.copy()
    #     ax.clear()
    #     maze_walls = np.ma.masked_where(I == 0, I)

    #     ax.imshow(
    #         c_copy[1:-1, 1:-1],
    #         origin="lower",
    #         extent=[0.0, Lx, 0.0, Ly],
    #         cmap="viridis",
    #         aspect="equal",
    #     )
    #     ax.imshow(
    #         maze_walls,
    #         origin="lower",
    #         extent=[0.0, Lx, 0.0, Ly],
    #         cmap="gray_r",
    #         alpha=0.15,
    #         interpolation="nearest",
    #         aspect="equal",
    #     )
    #     ax.set_xlabel("x")
    #     ax.set_ylabel("y")
    #     ax.set_title(f"Transient concentration, t = {n*dt:.2f}")
    #     plt.pause(0.01)
        

    # if n%50 == 0:
    #     print(n)



    #     # inlet & outlet coefficients of diffusion
    #     # c[1, 44] = c[1, 45] = c[1, 46] = c[1, 47] = c[1, 48] = 1
    #     # c[4,51] = c[5,51] =c[6,51] =c[7,51] =c[8,51] = 0

    
    if n == 1/dt or n == 5/dt or n == 15/dt or n == 50/dt: 
        plt.ion()
        c_copycat = c.copy()
        plt.clf()
        maze_walls = np.ma.masked_where(I == 0, I)
        plt.imshow(
            #no ghost cell plotting
            c_copycat[1:-1, 1:-1],
            origin="lower",
            extent=[0.0, Lx, 0.0, Ly],
            cmap="viridis",
        )
        plt.imshow(
            maze_walls,
            origin="lower",
            extent=[0.0, Lx, 0.0, Ly],
            cmap="gray_r",
            alpha=0.15,
            interpolation="nearest",
        )
        # plt.colorbar()
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(f"Transient concentration at T = {n*dt:.2f}")
        # plt.pause(ipause)
        plt.ioff()    
        plt.show() 
plt.ioff()
plt.close(fig)

end = time.time()
print('Implicit run finished.')
print(f"Simulation time: {end - start:.2f} seconds")


for j in range(1, Ny - 1):
    for i in range(1, Nx - 1):
        if (
        I[j, i] == 0
        and I[j, i + 1] == 0
        and I[j, i - 1] == 0
        and I[j + 1, i] == 0
        and I[j - 1, i] == 0
        ):
            dc_dx = (c[j+1, i+2] - c[j+1, i]) / (2*dx)
            dc_dy = (c[j+2, i+1] - c[j, i+1]) / (2*dy)
            diagnostic_quantity[j, i] = alpha_liquid * (np.sqrt(dc_dx**2 + dc_dy**2))
        else:
            diagnostic_quantity[j, i] = 0.0 
threshold = np.percentile(diagnostic_quantity, 90)
solution_mask = diagnostic_quantity > threshold
masked_solution = np.ma.masked_where(~solution_mask, diagnostic_quantity)
plt.figure()
maze_walls = np.ma.masked_where(I == 0, I)
plt.imshow(diagnostic_quantity, origin="lower", extent=[0.0, Lx, 0.0, Ly], cmap="inferno",
aspect="equal", vmin = 0, vmax = np.percentile(diagnostic_quantity, 95))
plt.imshow(maze_walls, origin="lower", extent=[0.0, Lx, 0.0, Ly], cmap="gray_r",
alpha=0.05, interpolation="nearest", aspect="equal")
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.show()

plt.figure(figsize=(6,5))
plt.imshow(c[1:Nx+1, 1:Ny+1].T, origin="lower", cmap="inferno")
plt.colorbar(label="c(x,y)")
plt.title("Final concentration field")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.show()

    
