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

def function_maker(route_of_administration, dose, t0, administration_time):
    if route_of_administration == 'inh':
        return inhalation_helper(dose, t0, administration_time), [t0, t0 + administration_time]
    elif route_of_administration == 'iv':
        return iv_infusion_helper(dose, t0, administration_time), [t0, t0 + administration_time]
    elif route_of_administration == 'oral':
        return oral_helper(dose, t0), [t0, t0 + 100*60]

def prepare_administration_functions(administration_regimen):
    functions = []
    bounds = []
    for administration_instance in administration_regimen:
        route_of_administration = administration_instance[0]
        dose = administration_instance[1]
        t0 = administration_instance[2]
        administration_time = administration_instance[3]
        function, bound = function_maker(route_of_administration, dose, t0, administration_time)
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

def administration_helper(administration_regimen):
    functions, bounds = prepare_administration_functions(administration_regimen)
    return whole_function_maker(functions, bounds), extract_times(administration_regimen)

class Administration:
    def __init__(self, route_of_administration, dose, t0=0, duration=0):
        self.dose = dose
        self.route_of_administration = route_of_administration
        self.t0 = t0
        self.duration = duration

    def __lt__(self, other):
        return self.t0 < other.t0

    def __str__(self):
        return f'ADMINISTRATION: ROA: {self.route_of_administration}, dose = {self.dose}, t0 = {self.t0}, duration = {self.duration}'

    def __add__(self, other):
        if type(other) == Administration:
            return Regimen(self, other)
        elif type(other) == Regimen:
            return other + self
        else:
            raise TypeError('Only Regimen and Administration are allowed')

    def repeat(self, n, t):
        output = Regimen(self)
        administration_instance = self
        for i in range(1, n):
            administration_instance = administration_instance.after(t)
            output += administration_instance
        return output

    def at(self, t0):
        dose = self.dose
        route_of_administration = self.route_of_administration
        duration = self.duration
        return Administration(route_of_administration, dose, t0, duration)

    def after(self, t):
        dose = self.dose
        route_of_administration = self.route_of_administration
        duration = self.duration
        t0 = self.t0
        return Administration(route_of_administration, dose, t0 + t, duration)       

class Regimen:
    def __init__(self, *administrations):
        self.list = [*administrations]
        self.list.sort()

    def __str__(self):
        output = 'REGIMEN:'
        for i, a in enumerate(self.list):
            output += f'\n   {i} ' + a.__str__()
        return output

    def __add__(self, other):
        if type(other) == Administration:
            self.list.append(other)
            return self
        elif type(other) == Regimen:
            self.list.extend(other.list)
            return self
        else:
            raise TypeError('Only Regimen and Administration are allowed')

    def repeat(self, n, t):
        output = self
        regimen_instance = self
        for i in range(1, n):
            regimen_instance = regimen_instance.after(t)
            output += regimen_instance
        return output           

    def after(self, t):
        output = Regimen()
        for i in self.list:
            output += i.after(t)
        return output        
        
    def at(self, t0):
        output = Regimen()
        origin = self.list[0].t0
        dt = t0 - origin
        for i in self.list:
            output += i.after(dt)
        return output

    def administration_function(self):
        return lambda t: 'da_dt'