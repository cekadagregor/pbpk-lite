import numpy as np
from scipy.optimize import minimize, Bounds

def interpolation(t_predicted, c_predicted, t_experimental):
    n = c_predicted.shape[0]
    result = np.empty((n, len(t_experimental)))
    for i in range(n):
        interpolation_instance = np.interp(t_experimental, t_predicted, c_predicted[i, -1, :])
        result[i, :] = interpolation_instance
    return result

def function_to_minimise(t_exp, cs_exp, x_to_t_and_c_function, x0):
    cs_exp = np.asarray(cs_exp)
    def inner(scaling):
        x = scaling * x0
        t, c = x_to_t_and_c_function(x)
        interpolated = interpolation(t, c, t_exp)
        mse = np.mean((interpolated - cs_exp) ** 2)
        return mse
    return inner

def fit_helper(t_exp, cs_exp, x_to_t_and_c_function, x0):
    f = function_to_minimise(t_exp, cs_exp, x_to_t_and_c_function, x0)
    output = minimize(f, np.ones(len(x0)), method='Nelder-Mead', bounds=Bounds(lb=1e-6, ub=np.inf))
    if output.success == False:
        return x0
    else:
        return output.x*x0

def fit_helper_oral(t_exp, cs_exp, x_to_t_and_c_function, x0):
    x0_cont = x0[:-1]
    t_max = 100*60
    results = []
    for n in range(2, 20):
        def f(scaling, n=n):
            x = scaling * x0_cont
            t, c = x_to_t_and_c_function(np.append(x, n))
            return np.mean((interpolation(t, c, t_exp) - cs_exp) ** 2)
        out = minimize(f, np.ones(len(x0_cont)), method='Nelder-Mead')
        results.append((out.fun, n, out.x * x0_cont))
    best = min(results, key=lambda r: r[0])
    return *best[2], best[1]
    


