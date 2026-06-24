"""
Main module for PBPK (Physiologically Based Pharmacokinetic) modeling.

This module provides the core model class for setting up and solving
PBPK models with customizable substance properties, patient parameters,
and elimination pathways.
"""

import numpy as np
from .distribution_params import V_function, Q_function, partition_model
from .graph import graph_whole_helper
from .solve import solver
from .ode_system import system_generator


class model:
    """
    Physiologically Based Pharmacokinetic (PBPK) model class.
    
    This class manages the setup and solution of PBPK models, handling
    substance properties, patient parameters, and drug elimination pathways.
    It integrates partition coefficients, blood flows, tissue volumes, and
    elimination kinetics to solve the system of ordinary differential equations.
    
    Attributes
    ----------
    partition_coefficients : ndarray
        Tissue-specific partition coefficients for the substance.
    blood_flows : ndarray
        Blood flow rates to each tissue.
    volumes : ndarray
        Volume of distribution for each tissue.
    elimination : callable
        Function that calculates elimination rate based on concentrations.
    odes : callable
        System of ordinary differential equations representing the model.
    t : ndarray
        Time points from the most recent solution.
    c : ndarray
        Concentrations in each tissue at corresponding time points.
    
    Examples
    --------
    >>> m = model()
    >>> m.set_substance(log_p=6.97, fu=0.0022448)
    >>> m.set_patient(bw=70)
    >>> doses = (10,)
    >>> times = (0, 24)
    >>> t, c = m.solve(doses, times)
    """
    def __init__(self):
        """
        Initialize the PBPK model with placeholder parameters.
        
        Sets up default substance properties (logp=6.97, fu=0.0022448),
        patient parameters (body weight=70 kg), and a default elimination
        pathway with hepatic and renal clearance.
        
        Notes
        -----
        The default parameters are placeholders and should be set using
        set_substance() and set_patient() methods before solving the model.
        Elimination is initialized with zero clearance values which can be
        updated using set_elimination().
        """
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
        Set the substance properties that are being investigated.

        Updates the partition coefficients based on the substance's
        physicochemical properties and regenerates the ODE system.

        Parameters
        ----------
        log_p : float
            Log of the octanol-water partition coefficient of the substance.
            Higher values indicate more lipophilic substances.
        fu : float
            Fraction unbound of the substance in blood (0 to 1).
            Represents the free unbound fraction available for distribution.
        
        Notes
        -----
        This method updates the partition_coefficients attribute and
        regenerates the ODE system with the new substance properties.
        The patient parameters must be set separately using set_patient().
        
        Examples
        --------
        >>> m = model()
        >>> m.set_substance(log_p=6.97, fu=0.0022448)
        """
    
        self.partition_coefficients = partition_model(log_p, fu)
        self.odes = system_generator(self.volumes, self.blood_flows, self.partition_coefficients, self.elimination)

    def set_patient(self, bw):
        """
        Set the patient's physiological distribution parameters.

        Updates blood flow rates and tissue volumes based on body weight,
        and regenerates the ODE system with the new patient parameters.

        Parameters
        ----------
        bw : float
            Patient's body weight in kilograms.
        
        Notes
        -----
        This method updates the blood_flows and volumes attributes by
        applying allometric scaling based on body weight. The ODE system
        is regenerated to reflect the new physiological parameters.
        Substance properties should be set separately using set_substance().
        
        Examples
        --------
        >>> m = model()
        >>> m.set_patient(bw=70)  # 70 kg patient
        """
        
        self.blood_flows = Q_function(bw)
        self.volumes = V_function(bw)
        self.odes = system_generator(self.volumes, self.blood_flows, self.partition_coefficients, self.elimination)
    
    def set_elimination(self, cl_l=0, cl_k=0):
        """
        Set the elimination parameters.

        Updates the elimination function to account for linear clearance
        from the liver and kidney. The elimination function calculates
        the outflow of substance from these tissues.

        Parameters
        ----------
        cl_l : float, optional
            Clearance from the liver based linearly on concentration
            (default is 0).
        cl_k : float, optional
            Clearance from the kidney based linearly on concentration
            (default is 0).
        
        Notes
        -----
        The elimination function modifies compartments 5 (liver) and 6 (kidney)
        with linear kinetics: outflow = clearance * concentration.
        After setting elimination parameters, the ODE system is regenerated.
        
        Examples
        --------
        >>> m = model()
        >>> m.set_elimination(cl_l=10, cl_k=5)
        """

        def elimination(c):
            outflow = np.zeros(16)
            outflow[5] += cl_k*c[5]
            outflow[6] += cl_l*c[6]
            return outflow

        self.elimination = elimination
        self.odes = system_generator(self.volumes, self.blood_flows, self.partition_coefficients, self.elimination)
    
    def solve(self, doses, times):
        """
        Set the administration regime and solve the system of ODEs.
        
        The `times` array must contain one more time point than the `doses`
        sequence. Each dose at index `i` is administered at `times[i]`, and
        the final entry in `times` is the last observation time.
        
        Parameters
        ----------
        doses : array-like
            A list of all doses.
        times : array-like
            A list of times at which doses are given, plus a final endpoint.
            This array must have length `len(doses) + 1`.
        
        Returns
        -------
        t : ndarray
            Times at which concentrations are calculated.
        c : ndarray
            Concentrations in each tissue at the given times.
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