from scipy.integrate import solve_ivp
import numpy as np

def solver(ode_system, times, volumes, n):
    t_full = np.array([])
    a_full = np.zeros((16*n,1))

    for i in range(len(times)):
        time_interval = times[i:i+2]
        a0 = a_full[:, -1]

        sol = solve_ivp(ode_system, time_interval, a0, method='Radau')

        t_partial = sol.t
        a_partial = sol.y

        t_full = np.concatenate((t_full, t_partial))
        a_full = np.concatenate((a_full, a_partial), axis=1)
        break

    
    a = np.delete(a_full.reshape(n, 16, -1), 0, 2)
    c_full = a/volumes.reshape(1,16,1)
    return t_full, c_full