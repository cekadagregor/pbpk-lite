"""
Main module for PBPK (Physiologically Based Pharmacokinetic) modeling.

This module provides the core model class for setting up and solving
PBPK models with customizable substance and metabolite properties,
patient parameters, and elimination pathways.
"""

from .distribution_params import V_function, distribution_helper
from .graph import graph_whole_helper, graph_venous_helper, graph_compartments_helper
from .solve import solver
from .ode_system import system_generator
from .administration import administration_helper
from .metabolism_and_excretion import metabolism_helper
from .optimisation import fit_helper, fit_helper_oral
import numpy as np

class model:
    """
    Physiologically Based Pharmacokinetic (PBPK) model class.
    
    This class manages the setup and solution of PBPK models, handling
    substance properties, patient parameters, and drug elimination pathways.
    It integrates partition coefficients, blood flows, tissue volumes, and
    elimination kinetics to solve the system of ordinary differential equations.
    The implementation uses nanograms (ng) for doses and amounts, milliliters
    (mL) for volumes, milliliters per minute (mL/min) for blood flows and
    clearances, nanograms per milliliter (ng/mL) for concentrations, and
    minutes for time.
    
    Attributes
    ----------
    partition_coefficients : list of ndarray
        Tissue-specific partition coefficients for each substance.
    blood_flows : ndarray
        Blood flow rates to each tissue.
    volumes : ndarray
        Volume of distribution for each tissue.
    eliminations : list of callable
        Functions that calculate elimination rates for each substance.
    gravimetric_coefficients : list of float
        Molecular-mass ratios used to transfer eliminated parent substance to
        its metabolites.
    odes : callable
        System of ordinary differential equations representing the model.
    t : ndarray
        Time points from the most recent solution.
    c : ndarray
        Concentrations in each compartment and substance at corresponding time
        points. Rows are grouped in blocks of 16 compartments per substance.
    
    Examples
    --------
    >>> m = model()
    >>> m.set_substance_and_metabolites(
    ...     (6.97, 5.33, 5.24), (0.0022448, 0.01209, 0.101),
    ...     (314.469, 330.468, 344.451))
    >>> m.set_patient(bw=70)
    >>> m.set_administration([('iv', 0.053, 0, 2)])
    >>> m.set_metabolism_and_excretion(
    ...     ((0, 1), (1, 2), (2, -1), (2, -1)),
    ...     ((6, 247), (6, 998), (6, 8.94), (5, 1.48)))
    >>> t, c = m.simulate(additional_time=24 * 60)
    """
    default_oral_params=(20, 450, 2)
    def __init__(self):
        """
        Initialize the PBPK model with placeholder parameters.
        
        Sets up default distribution parameters for a 70 kg patient.
        
        Notes
        -----
        The default parameters are example values and can be replaced using
        set_substance_and_metabolites(), set_patient(),
        set_administration(), and set_metabolism_and_excretion().
        """
        # placeholder parameters
        logps = (6.97, 5.33, 5.24)
        self.n = len(logps)
        fus = (0.0022448, 0.01209, 0.101)  
        
        bw = 70
        self.body_weight = bw
        mms = (314.469, 330.468, 344.451)
        con = ((0, 1), (1, 2), (2, -1), (2, -1))
        cls = ((6, 247), (6, 998), (6, 8.94), (5, 1.48))
        self.volumes = V_function(bw)

    def set_substance_and_metabolites(self, logps, fus, mms):
        """
        Set the substance properties that are being investigated.

        Updates partition coefficients and metabolite mass ratios based on the
        physicochemical properties and molecular masses.

        Parameters
        ----------
        logps : sequence of float
            Log octanol-water partition coefficients for the substance and
            metabolites.
        fus : sequence of float
            Fractions unbound in blood, with one value for each substance.
        mms : sequence of float
            Molecular masses, with one value for each substance. The sequence
            order defines the parent substance followed by its metabolites.
        
        Notes
        -----
        This method updates the substance and distribution parameters used by
        the model.
        The three sequences must have the same nonzero length.
        
        Examples
        --------
        >>> m = model()
        >>> m.set_substance_and_metabolites(
        ...     (6.97, 5.33, 5.24), (0.0022448, 0.01209, 0.101),
        ...     (314.469, 330.468, 344.451))
        """

        self.molar_masses = mms
        self.n = len(logps)
        self.substance_parameters = logps, fus
        self.distribution = distribution_helper(self.body_weight, logps, fus)

    def set_patient(self, bw):
        """
        Set the patient's physiological distribution parameters.

        Updates blood flow rates and tissue volumes based on body weight.

        Parameters
        ----------
        bw : float
            Patient's body weight in kilograms.
        
        Notes
        -----
        This method updates the volumes and distribution parameters by
        applying allometric scaling based on body weight.
        Substance and metabolite properties can be set separately using
        set_substance_and_metabolites().
        
        Examples
        --------
        >>> m = model()
        >>> m.set_patient(bw=70)  # 70 kg patient
        """

        self.body_weight = bw
        self.volumes = V_function(bw)
        self.distribution = distribution_helper(bw, *self.substance_parameters)

    def set_administration(self, regimen, oral_params=default_oral_params):
        """
        Set the administration regimen.

        Parameters
        ----------
        regimen : sequence of sequence
            Administration entries in the form
            ``(route, dose, start_time, administration_time)``. Times and
            doses are expressed in minutes and nanograms, respectively.
            Supported routes are ``'iv'``, ``'inh'``, and ``'oral'``.
        """
        self.latest_regimen = regimen
        self.administration, self.times = administration_helper(regimen, oral_params)
                    
    def set_metabolism_and_excretion(self, connections, clearances):
        """
        Set the metabolism and excretion parameters.

        Updates the metabolism and excretion function to account for linear
        clearance from the liver and kidney.

        Parameters
        ----------
        connections : sequence of sequence of int
            Parent-to-metabolite connections. Each pair contains the source
            substance index and the destination index; ``-1`` denotes
            excretion from the system.
        clearances : sequence of sequence of float
            Clearance entries in the form ``(compartment, clearance)``.
        
        Notes
        -----
        Clearance is modeled with linear kinetics:
        ``outflow = clearance * concentration``.
        
        Examples
        --------
        >>> m = model()
        >>> m.set_metabolism_and_excretion(
        ...     ((0, 1), (1, 2), (2, -1)),
        ...     ((6, 10), (6, 10), (5, 5)))
        """

        self.connections = np.array(connections)
        self.clearances = np.array(clearances)
        self.metabolism_and_excretion = metabolism_helper(self.molar_masses, connections, clearances)

    def simulate(self, additional_time=60):
        """
        Simulate the configured system of ordinary differential equations.

        Call set_administration() before simulation. The simulation continues
        for ``additional_time`` minutes after the final configured
        administration time.
        
        Parameters
        ----------
        additional_time : float, optional
            Number of minutes to simulate after the final administration time.
            The default is 60.
        
        Returns
        -------
        t : ndarray
            Times at which concentrations are calculated.
        c : ndarray
            Concentrations in each compartment at the given times. Rows are
            grouped in blocks of 16 compartments per substance.
        
        Examples
        --------
        >>> m = model()
        >>> m.set_substance_and_metabolites(
        ...     (6.97, 5.33, 5.24), (0.0022448, 0.01209, 0.101),
        ...     (314.469, 330.468, 344.451))
        >>> m.set_patient(bw=70)
        >>> m.set_administration([('iv', 0.053, 0, 2)])
        >>> m.set_metabolism_and_excretion(
        ...     ((0, 1), (1, 2), (2, -1)),
        ...     ((6, 247), (6, 998), (5, 1.48)))
        >>> t, c = m.simulate(additional_time=24 * 60)
        """

        odes = system_generator(self.volumes, self.n, self.administration, self.distribution, self.metabolism_and_excretion)
        times = np.append(self.times, additional_time + self.times[-1])
        self.t, self.c = solver(odes, times, self.volumes, self.n)
        return self.t, self.c

    def fit_clearances(self, t, c, report=True):
        '''best if t, c is iv bolus experimental data'''
        def add_positions(cl_only):
            clearances = self.clearances
            clearances[:, 1] = cl_only
            return clearances

        def f(clearances_only):
            self.set_metabolism_and_excretion(self.connections, add_positions(clearances_only))
            return self.simulate(t[-1])

        new_clearances = fit_helper(t, c, f, self.clearances[:, 1])
        if report == True:
            print(f'the new clearances are {new_clearances}')
        self.metabolism_and_excretion = metabolism_helper(self.molar_masses, self.connections, add_positions(new_clearances))

    def fit_oral(self, t, c, report=True):
        '''best if you first fit iv bolus'''
        def f(oral_params):
            self.set_administration(self.latest_regimen, oral_params)
            return self.simulate(t[-1])

        new_oral_params = fit_helper_oral(t, c, f, self.default_oral_params)
        if report == True:
            print(f'the new clearances are {new_oral_params}')
        self.default_oral_params = new_oral_params
        
    def graph_whole(self, name, time_unit='min'):
        """
        Graph the concentrations in each tissue and save the figure.
        
        Parameters
        ----------
        name : str
            File path to save the figure.
        time_unit : {'min', 'hours', 'days'}, optional
            Unit used for the x-axis label and values. The default is 'min'.
        """
        if self.n == 1:
            graph_whole_helper(self.t, self.c, name, time_unit=time_unit)
        else:
            for i in range(self.n):
                new_name = name + str(i)
                graph_whole_helper(self.t, self.c, new_name, time_unit=time_unit, substance_index=i)         

    def graph_venous(self, name, limit_of_detection=None, log=True, time_unit='min'):
        """
        Graph the concentrations in the venous blood compartment and save the figure.

        Parameters
        ----------
        name : str
            File path to save the figure.
        limit_of_detection : float, optional
            Limit of detection to be indicated on the graph.
        log : bool, optional
            Whether to use a logarithmic y-axis. The default is True.
        time_unit : {'min', 'hours', 'days'}, optional
            Unit used for the x-axis label and values. The default is 'min'.
        """
        graph_venous_helper(self.t, self.c, name, limit_of_detection, log, time_unit=time_unit)

    def graph_compartments(self, compartments, name, time_unit='min'):
        """
        Graph the concentrations in selected compartments and save the figure.

        Parameters
        ----------
        compartments : list of str or int
            List of compartment names or indices to be graphed.
        name : str
            File path to save the figure.
        time_unit : {'min', 'hours', 'days'}, optional
            Unit used for the x-axis label and values. The default is 'min'.
        """
        graph_compartments_helper(self.t, self.c, compartments, name, time_unit=time_unit)