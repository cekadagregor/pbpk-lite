import numpy as np

def system_generator(V, n, administration, distribution, metabolism_and_excretion):
    def inner(t, A):
        dA = np.zeros((n, 16))
        C = A.reshape((n, 16))/V

        dA[0, :] += administration(t)
        dA += distribution(C)
        dA += metabolism_and_excretion(C)
        return dA.flatten()
    return inner