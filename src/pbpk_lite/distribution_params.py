import numpy as np

V_a = np.array((143, 124, 20.7, 23.6, 3.8, 4.4, 24.1, 16.7, 429, 1.2, 111, 2.7, 2.2, 0.51, 25.7, 51.4)) # mL/kg
def V_function(bw): return V_a * bw

Q_a = np.array((3.7, 3.6, 10, 13, 2.14, 15.7, 21, 71, 10.7, 1.9, 4.3, 1.1, 0.56, 0.04)) # mL/min/kg
def Q_function(bw): return Q_a * (bw**0.75)

# phospholipid, neutral lipid, water
Tissue_composition = np.array(((0.002, 0.79, 0.18),
(0.0005, 0.074, 0.439),
(0.0565, 0.051, 0.77),
(0.0163, 0.0487, 0.718),
(0.0166, 0.0115, 0.758),
(0.0162, 0.0207, 0.783),
(0.0252, 0.0348, 0.751),
(0.009, 0.003, 0.811),
(0.0072, 0.0238, 0.76),
(0.0188, 0.0723, 0.66),
(0.0111, 0.0284, 0.718),
(0.0198, 0.0201, 0.788),
(0.0182, 0.0338, 0.784),
(0, 0, 0.859),
(0.0033, 0.0022, 0.651),
(0.00225, 0.0035, 0.945))) # L/kg

f_pl = Tissue_composition[:, 0]
f_nl = Tissue_composition[:, 1]
f_w = Tissue_composition[:, 2]

f_pl_p = 0.067
f_nl_p = 0.0012
f_w_p = 0.939

def partition_model(logp, fu):
    fut = 1/(1 + 0.5*((1-fu)/fu))
    p_ow = 10**logp
    # from https://doi.org/10.1002/jps.10005
    log_d_ow = 1.115*logp - 1.35
    d_ow = 10**log_d_ow


    k_p = (p_ow*(f_nl + 0.3*f_pl) + (f_w + 0.7*f_pl)/fut)/(p_ow*(f_nl_p + 0.3*f_pl_p) + (f_w_p + 0.7*f_pl_p)/fu)
    k_p[0] = (d_ow*(f_nl[0] + 0.3*f_pl[0]) + f_w[0] + 0.7*f_pl[0])/(d_ow*(f_nl_p + 0.3*f_pl_p) + (f_w_p + 0.7*f_pl_p)/fu)
    return k_p