import numpy as np

def metabolism_instance(molar_masses, connection, clearance_array, C, n):
    dA = np.zeros((n, 16))
    flow_from = connection[0]
    flow_to = connection[1]
    flow = clearance_array * C[flow_from, :]
    dA[flow_from, :] -= flow
    if flow_to != -1: #-1 index means excretion
        dA[flow_to, :] += flow * molar_masses[flow_to]/molar_masses[flow_from]
    return dA

def preprocesing(clearances):
    clearances = np.array(clearances)
    n, l = clearances.shape
    if l == 16:
        return clearances
    clearances_new = np.zeros((n, 16))
    for connection_index in range(n):
        clearance_pair = clearances[connection_index, :]
        organ_index = int(clearance_pair[0])
        clearance = clearance_pair[1]
        clearances_new[connection_index, organ_index] = clearance
    return clearances_new

def metabolism_helper(molar_masses, connections, clearances):
    n = len(molar_masses)
    clearances = preprocesing(clearances)
    def inner(C):
        dA = np.zeros((n, 16))
        for i, connection in enumerate(connections):
            clearance_array = clearances[i, :]
            dA += metabolism_instance(molar_masses, connection, clearance_array, C, n)
        return dA
    return inner
