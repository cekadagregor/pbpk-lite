import numpy as np
import functools

max_whole_gut_transit_time = 100*60

def curry_t(mechanism):
    @functools.wraps(mechanism)
    def inner(**kwargs):
        return lambda t: mechanism(t, **kwargs)
    return inner

@curry_t
def inhalation_mechanism(t, dose, duration):
    flow_vector = np.zeros(16)
    flow = dose/duration
    flow_vector[7] += flow
    return flow_vector

@curry_t
def oral_mechanism(t, dose, duration = max_whole_gut_transit_time, t_lag=20, scaling_factor=450, shape_parameter=2):
    flow_vector = np.zeros(16)
    weibull_distribution_derivative = -(dose*shape_parameter/scaling_factor) * ((t-t_lag)/scaling_factor)**(shape_parameter-1) * np.exp(-((t-t_lag)/scaling_factor)**shape_parameter)
    if weibull_distribution_derivative < 0:
        flow_vector[6] += -weibull_distribution_derivative
    return flow_vector

@curry_t
def iv_infusion_mechanism(t, dose, duration):
    flow_vector = np.zeros(16)
    flow = dose/duration
    flow_vector[15] += flow
    return flow_vector

def helper_constructor(route_of_administration):
    if route_of_administration == 'inh':
        return inhalation_mechanism
    elif route_of_administration == 'oral':
        return oral_mechanism
    elif route_of_administration == 'iv':
        return iv_infusion_mechanism

class Administration:
    def __init__(self, route_of_administration, dose, t0=0, duration=0):
        self.dose = dose
        self.route_of_administration = route_of_administration
        self.t0 = t0
        if route_of_administration == 'oral' and duration == 0:
            duration = max_whole_gut_transit_time
        self.duration = duration
        self.helper = helper_constructor(route_of_administration)

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

    def administration_function(self):
        f_of_t = self.helper(dose = self.dose, duration = self.duration)
        def inner(t):
            output = np.zeros(16)
            if self.t0 <= t < self.t0 + self.duration:
                output += f_of_t(t)
            return output
        return inner


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
        def inner(t):
            def each(index):
                administration_instance = self.list[index]
                f = administration_instance.administration_function()
                if index == -1:
                    return f(t)
                else:
                    return f(t) + each(index - 1)
            return each(0)
        return inner