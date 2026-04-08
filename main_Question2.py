import numpy as np
import matplotlib.pyplot as plt
import time
def solver_loop(N, animate = True):
    Lx = Ly = 1.0
    Nx = Ny = N
    dx = Lx / Nx
    dy = Ly / Ny

    alpha = 1.0e-3
    U0 = 1.0
    T_left = 1.0
    T_right = 0.0

    nstep = 20000
    if N < 50:
        dt = 0.02
    elif N < 70:
        dt = 0.010
    else: 
        dt = 0.0075

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
    if animate == True:
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

    u_plus = np.maximum(u, 0.0)
    u_minus = np.minimum(u, 0.0)
    v_plus = np.maximum(v, 0.0)
    v_minus = np.minimum(v, 0.0)

    start = time.time()
    for istep in range(nstep):
        time_step += dt

        T_new = T.copy()

        T[:, 0] = T_left
        T[:, Nx+1] = T_right
        T[0, :] = T[1, :]
        T[Ny+1, :] = T[Ny, :]

        #comment the for loops for this solver
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
                dif = alpha * (((T[j,i+1] - 2*T[j,i] + T[j,i-1])/dx**2) + ((T[j+1,i] - 2*T[j,i] + T[j-1,i])/dy**2))

                # Diffusion - Advection

                T_new[j,i] = T[j,i] - dt*(adv_x +adv_y) + dt*dif
            
        
        T_new[:, 0] = T_left
        T_new[:, Nx+1] = T_right
        T_new[0, :] = T_new[1, :]
        T_new[Ny+1, :] = T_new[Ny, :]

        delta = np.max(np.abs(T_new-T))
        # Print nstep and time
        if istep % 1000 == 0:
            print(f'Time ={time_step:5.2f} seconds | Iteration={istep:3d} | Delta = {delta}')
        

        # Steady State

        if delta < 1e-6:
            if animate == True:
                img.set_data(T[1:-1, 1:-1])
                plt.pause(0.001)
            print('Steady state reached')
            print(f'Iteration: {istep}')
            print(f'Time: {time_step:.2f} seconds')
            print(f'Delta:{delta:.3e}')
            T = T_new.copy()
            break

        T = T_new.copy()

        if animate and istep % 8000 == 0:   # update every 20 steps for speed
            img.set_data(T[1:-1, 1:-1])
            plt.pause(0.001)

    end = time.time()
    print('Explicit run finished.')
    print(f"Simulation time: {end - start:.2f} seconds")

    mid_x = Nx//2
    mid_y = Ny//2
    T_noGhost = T[1:-1, 1:-1]

    #Centerline Profile X-axis
    centerline_x= T_noGhost[:,mid_x]

    #Centerline Profile Y-axis
    centerline_y = T_noGhost[mid_y,:]
    if animate:
        plt.ioff()
        plt.show()
    return y, x, centerline_y, centerline_x



N_values = [20, 40, 60, 80]

def run_solver_set(solver_func, N_values):
    results = {}
    for N in N_values:
        y, x, centerline_y, centerline_x = solver_func(N, animate=False)
        results[N] = {
            "xgrid": x,
            "ygrid": y,
            "centerline_x": centerline_x,
            "centerline_y": centerline_y,
        }
    return results 

results_loop = run_solver_set(solver_loop, N_values)


plt.figure()
for N in N_values:
    plt.plot(results_loop[N]["ygrid"], results_loop[N]["centerline_x"], label=f"loop, N={N}")
 
plt.xlabel("y")
plt.ylabel("T")
plt.title("Centerline profile at x = 0.5")
plt.legend()
plt.grid()
plt.show()

plt.figure()
for N in N_values:
    plt.plot(results_loop[N]["xgrid"], results_loop[N]["centerline_y"], label=f"loop, N={N}")
plt.xlabel("x")
plt.ylabel("T")
plt.title("Centerline profile at y = 0.5")
plt.legend()
plt.grid()
plt.show()

plt.ioff()
plt.show()



