import numpy as np
from numba import njit

# 0 adipose, 1 bone, 2 brain, 3 gut, 4 heart, 5 kidney, 6 liver, 7 lung, 8 muscle, 
# 9 pancreas, 10 skin, 11 spleen, 12 stomach, 13 testes, 14 arterial blood, 15 venous blood

arterial_input = np.ones(14)
arterial_input[7] = 0

liver_input = np.zeros(14)
liver_input[11] = 1
liver_input[9] = 1
liver_input[12] = 1
liver_input[3] = 1

venous_input = np.ones(14)
venous_input[7] = 0
venous_input[11] = 0
venous_input[9] = 0
venous_input[12] = 0
venous_input[3] = 0

no_lung = np.ones(14)
no_lung[7] = 0

iv_input = 0
ia_input = 0

# @njit
def system_generator(V, Q, K, elimination):
    def inner(t, A):
        C = A/V
        dA = np.zeros(16)
        dA[7] = Q[7]*(C[15] - C[7]/K[7])
        dA[:14] += no_lung*Q*(C[14] - C[:14]/K[:14])
        dA[6] += sum(liver_input*Q*C[:14]/K[:14])
        dA -= elimination(C)
        dA[15] = - Q[7]*C[15] + sum(venous_input*Q*C[:14]/K[:14]) + iv_input
        dA[14] = Q[7]*C[7]/K[7] - (sum(Q) - Q[7])*C[14] + ia_input
        return dA
    return inner