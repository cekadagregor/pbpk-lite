import numpy as np

def clearance_helper(cl_l, cl_k):
    def inner(c):
        outflow = np.zeros(16)
        outflow[5] += cl_k * c[5]
        outflow[6] += cl_l * c[6]
        return outflow
    return inner