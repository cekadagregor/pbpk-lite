import numpy as np
from scipy.optimize import minimize

def interpolation(t_predicted, c_predicted, t_experimental):
    n = c_predicted.shape[0]
    result = np.empty((n, len(t_experimental)))
    for i in range(n):
        interpolation_instance = np.interp(t_experimental, t_predicted, c_predicted[i, -1, :])
        result[i, :] = interpolation_instance
    return result

def mean_square_error(c_predicted, c_experimental):
    square_error = (c_predicted - c_experimental)**2
    return np.mean(square_error)

def function_to_minimise(t_exp, c_exp, x_to_t_and_c_function):
    def inner(x):
        t, c = x_to_t_and_c_function(x)
        interpolated = interpolation(t, c, t_exp)
        mse = mean_square_error(interpolated, c_exp)
        return mse
    return function_to_minimise

def fit_helper(t_exp, c_exp, x_to_t_and_c_function, x0):
    f = function_to_minimise(t_exp, c_exp, x_to_t_and_c_function)
    output = minimize(f, x0)
    if output.succes == False:
        return x0
    else:
        return output.x


