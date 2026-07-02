from scipy.integrate import solve_ivp
import numpy as np

def solver(ode_system, doses, times, volumes, route_of_administration):
    if len(times) != len(doses) + 1:
        raise ValueError("times must contain exactly len(doses) + 1 entries")

    if np.any(np.diff(times) <= 0):
        raise ValueError("times must be strictly increasing")

    if route_of_administration not in {'iv', 'ia', 'inh'}:
        raise ValueError("route_of_administration must be one of 'iv', 'ia', or 'inh'")

    t_full = np.array([])
    a_full = np.array([[] for i in range(16)])

    for i in range(len(doses)):
        dose = doses[i]
        time_interval = times[i:i+2]
        if i == 0:
            a0 = np.zeros(16)
        else:
            a0 = a_full[:, -1]

        if route_of_administration == 'iv':
            index = 15
        elif route_of_administration == 'ia':
            index = 14
        elif route_of_administration == 'inh':
            index = 7
        a0[index] += dose

        sol = solve_ivp(ode_system, time_interval, a0, method='Radau')
        t_partial = sol.t
        a_partial = sol.y

        t_full = np.concatenate((t_full, t_partial))
        a_full = np.concatenate((a_full, a_partial), axis=1)

    volumes = volumes[:, np.newaxis]
    c_full = a_full/volumes
    return t_full, c_full