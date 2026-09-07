import numpy as np

def inhalation_helper(dose, t0, administration_time):
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