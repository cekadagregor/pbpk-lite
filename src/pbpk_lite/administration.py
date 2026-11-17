import numpy as np

def inhalation_helper(dose, t0, administration_time=10):
    def inner(t):
        output = np.zeros(16)
        t -= t0
        da_dt = dose/administration_time
        if t <= administration_time and t >= 0:
            output[7] += da_dt
        return output
    return inner

def placeholder_helper():
    def inner(t):
        return np.zeros(16)
    return inner

def oral_helper(dose, t0, t_lag=20, scaling_factor=450, shape_parameter=2):
    def inner(t):
        output = np.zeros(16)
        t_max = 100*60 #100 h max whole gut transit time
        t -= t0
        if t <= t_max:
            weibull_distribution_derivative = -(dose*shape_parameter/scaling_factor) * ((t-t_lag)/scaling_factor)**(shape_parameter-1) * np.exp(-((t-t_lag)/scaling_factor)**shape_parameter)
            if weibull_distribution_derivative < 0:
                output[6] += -weibull_distribution_derivative
        return output
    return inner

def iv_infusion_helper(dose, t0, administration_time):
    def inner(t):
        output = np.zeros(16)
        t -= t0
        da_dt = dose/administration_time
        if t <= administration_time:
            output[15] += da_dt
        return output
    return inner

def function_maker(route_of_administration, dose, t0, administration_time, oral_params):
    if route_of_administration == 'inh':
        return inhalation_helper(dose, t0, administration_time), [t0, t0 + administration_time]
    elif route_of_administration == 'iv':
        return iv_infusion_helper(dose, t0, administration_time), [t0, t0 + administration_time]
    elif route_of_administration == 'oral':
        return oral_helper(dose, t0, *oral_params), [t0, t0 + 100*60]

def prepare_administration_functions(administration_regimen, oral_params):
    functions = []
    bounds = []
    for administration_instance in administration_regimen:
        route_of_administration = administration_instance[0]
        dose = administration_instance[1]
        t0 = administration_instance[2]
        administration_time = administration_instance[3]
        function, bound = function_maker(route_of_administration, dose, t0, administration_time, oral_params)
        functions.append(function)
        bounds.append(bound)
    return functions, bounds

def extract_times(administration_regimen):
    administration_regimen = np.array(administration_regimen)
    times = administration_regimen[:, 2]
    times = times.astype(float)
    return times

def whole_function_maker(functions, bounds):
    def inner(t):
        output = np.zeros(16)
        for i, bound in enumerate(bounds):
            if bound[0] <= t <= bound[1]:
                output += functions[i](t)
        return output
    return inner

def administration_helper(administration_regimen, oral_params):
    functions, bounds = prepare_administration_functions(administration_regimen, oral_params)
    return whole_function_maker(functions, bounds), extract_times(administration_regimen)
