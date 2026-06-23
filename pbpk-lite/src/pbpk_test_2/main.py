import numpy as np
from .distribution_params import V_function, Q_function, partition_model
from .graph import graph_whole_helper
from .solve import solver
from .ode_system import system_generator

class model:
    def __init__(self):
        # placeholder parameters
        logp = 6.97
        fu = 0.0022448
        bw = 70
        self.partition_coefficients = partition_model(logp, fu)
        self.blood_flows = Q_function(bw)
        self.volumes = V_function(bw)

        cl_k = 0
        cl_l = 10
        def elimination(c):
            outflow = np.zeros(16)
            outflow[5] += cl_k*c[5]
            outflow[6] += cl_l*c[6]
            return outflow

        self.elimination = elimination

    def set_substance(self, log_p, fu):
        """
        Set the substance that is being investigated.

        Parameters
        ----------
        log_p : float
            Log of the octanol water partition of the substance.
        fu : float
            Free unbound value of the substance.
        """
    
        self.partition_coefficients = partition_model(log_p, fu)
        self.odes = system_generator(self.volumes, self.blood_flows, self.partition_coefficients, self.elimination)

    def set_patient(self, bw):
        """
        Set the patient's distribution parameters.

        Parameters
        ----------
        bw : float
            Patient's body weight.
        """
        
        self.blood_flows = Q_function(bw)
        self.volumes = V_function(bw)
        self.odes = system_generator(self.volumes, self.blood_flows, self.partition_coefficients, self.elimination)
    
    def set_elimination(self, /, cl_l = 0, cl_k= 0):
        '''
        Set the elimination parameters.

        Parameters
        ----------
        cl_l : float
            Clearance from the liver based linearly on concentration.
        cl_k : float
            Clearance from the kidney based linearly on concentration.
        '''

        def elimination(c):
            outflow = np.zeros(16)
            outflow[5] += cl_k*c[5]
            outflow[6] += cl_l*c[6]
            return outflow

        self.elimination = elimination
        self.odes = system_generator(self.volumes, self.blood_flows, self.partition_coefficients, self.elimination)
    
    def solve(self, doses, times):
        """
        Set the administration regime and solve the system of odes.
        
        Each dose at index i is administered at time at index i.
        
        Parameters
        ----------
        doses : array-like
            A list of all doses.
        times : array-like
            A list of times at which doses ware given.
        
        Returns
        -------
        t : ndarray
            Times at which there is some c.
        c : ndarray
            Concentrations in each tissue at given t.
        """

        self.t, self.c = solver(self.odes, doses, times, self.volumes)
        return self.t, self.c
        
    def graph_whole(self, name):
        """
        Graph the concentrations in each tissue and save the figure.
        
        Parameters
        ----------
        name : str
            File path to save the figure.
        """
        graph_whole_helper(self.t, self.c, name)