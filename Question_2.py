import numpy as np
import matplotlib.pyplot as plt
import time

Lx = Ly = 1.0
Nx = Ny = 20
dx = Lx / Nx
dy = Ly / Ny

alpha = 1.0e-3
U0 = 1.0
T_left = 1.0
T_right = 0.0

nstep = 50000
dt = 0.02

x = np.linspace(0.5 * dx, Lx- 0.5 * dx, Nx)
y = np.linspace(0.5 * dy, Ly- 0.5 * dy, Ny)
X, Y = np.meshgrid(x, y)

u =-U0 * np.sin(np.pi * X) * np.cos(np.pi * Y)
v = U0 * np.cos(np.pi * X) * np.sin(np.pi * Y)

T = np.full((Ny + 2, Nx + 2), T_right)

time_step = 0.0

print('Courant Number: {}'.format(U0*dt/dx))
print('Cell-based Peclet number: {}'.format(U0*dx/alpha))
print('"Diffusive" Courant number: {}'.format(dt*alpha/dx**2))

plt.ion()
fig, ax = plt.subplots()

img = ax.imshow(
    T[1:-1, 1:-1],
    origin="lower",
    extent=[0, Lx, 0, Ly],
    cmap="inferno",
    aspect="equal",
    vmax = 1,
    vmin = 0
)

quiv = ax.quiver(X, Y, u, v, color="white")
plt.colorbar(img, ax=ax)
plt.pause(0.01)

start = time.time()
for istep in range(nstep):
    time_step += dt

    T_new = T.copy()

    T[:, 0] = T_left
    T[:, Nx+1] = T_right
    T[0, :] = T[1, :]
    T[Ny+1, :] = T[Ny, :]

    for j in range(1, Ny + 1):
        for i in range(1, Nx+1):

            # +ve and -ve velocities
            u_plus = max(u[j-1,i-1], 0)
            u_minus = min(u[j-1,i-1], 0)
            v_plus = max(v[j-1,i-1], 0)
            v_minus = min(v[j-1,i-1], 0)

            # Advection x 
            adv_x = (u_plus * (T[j,i] - T[j,i-1]) + u_minus * (T[j,i+1] - T[j,i]))/dx

            # Advection y

            adv_y = (v_plus * (T[j,i] - T[j-1,i]) + v_minus * (T[j+1,i] - T[j,i]))/dy

            #Central Diffusion 
            dif = alpha * (((T[j,i+1] - 2*T[j,i] + T[j,i-1])/dx**2)+ ((T[j+1,i] - 2*T[j,i] + T[j-1,i])/dy**2))

            # Diffusion - Advection

            T_new[j,i] = T[j,i] - dt*(adv_x +adv_y) + dt*dif

    T_new[:, 0] = T_left
    T_new[:, Nx+1] = T_right
    T_new[0, :] = T_new[1, :]
    T_new[Ny+1, :] = T_new[Ny, :]

    # Print nstep and time
    if istep % 1000 == 0:
        print(f'Time ={time_step:5.2f} seconds | Iteration={istep:3d}')
    

    # Steady State
    delta = np.max(np.abs(T_new-T))
    if delta < 1e-6:
        img.set_data(T[1:-1, 1:-1])
        plt.pause(0.001)
        print('Steady state reached')
        print(f'Iteration: {istep}')
        print(f'Time: {time_step:.2f} seconds')
        print(f'Delta:{delta:.3e}')
        T = T_new.copy()
        break

    T = T_new.copy()

    if istep % 2000 == 0:   # update every 20 steps for speed
        img.set_data(T[1:-1, 1:-1])
        plt.pause(0.001)

end = time.time()
print('Explicit run finished.')
print(f"Simulation time: {end - start:.2f} seconds")


mid_x = Nx//2
mid_y = Ny//2
T_noGhost = T[1:-1, 1:-1]

#Centerline Profile X-axis
centerline_x = T_noGhost[:,mid_x]

#Centerline Profile Y-axis
centerline_y = T_noGhost[mid_y,:]


plt.ioff()
# plt.figure()
# plt.imshow(
# T[1:-1, 1:-1],
# origin="lower",
# extent=[0.0, Lx, 0.0, Ly],
# cmap="inferno",
# aspect="equal"
# )
# plt.colorbar()
# plt.quiver(X, Y, u, v, color="white")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.show()

plt.figure()
plt.plot(y, centerline_x)
plt.xlabel("y")
plt.ylabel("T_x")
plt.show(block=False)

plt.figure()
plt.plot(x, centerline_y)
plt.xlabel("x")
plt.ylabel("T_y")
plt.show()