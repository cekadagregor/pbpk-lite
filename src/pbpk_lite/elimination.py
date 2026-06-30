import numpy as np
from numba import njit

@njit
def elimin_func(cl_l, cl_k, c):
    outflow = np.zeros(16)
    outflow[5] += cl_k * c[5]
    outflow[6] += cl_l * c[6]
    return outflow  

def clearance_helper(cl_l, cl_k):
    def inner(c):
        return elimin_func(cl_l, cl_k, c)
    return inner