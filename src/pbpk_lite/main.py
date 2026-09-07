"""
Main module for PBPK (Physiologically Based Pharmacokinetic) modeling.

This module provides the core model class for setting up and solving
PBPK models with customizable substance and metabolite properties,
patient parameters, and elimination pathways.
"""

from .distribution_params import V_function, Q_function, partitions, gravimetric_model
from .graph import graph_whole_helper, graph_venous_helper, graph_compartments_helper
from .solve import solver
from .ode_system import system_generator
from .elimination import clearance_helper
from .administration import inhalation_helper, placeholder_helper

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
    >>> m.set_elimination((247, 998, 8.94), (0, 0, 1.48))
    >>> doses = [0.053]
    >>> times = [0, 24*60]
    >>> t, c = m.simulate(doses, times)
    """
    def __init__(self):
        """
        Initialize the PBPK model with placeholder parameters.
        
        Sets up default properties for one substance and two metabolites,
        patient parameters (body weight=70 kg), and linear hepatic and renal
        clearance for each substance.
        
        Notes
        -----
        The default parameters are example values and can be replaced using
        set_substance_and_metabolites(), set_patient(), and set_elimination().
        """
        # placeholder parameters
        logps = (6.97, 5.33, 5.24)
        fus = (0.0022448, 0.01209, 0.101)  
        mms = (314.469, 330.468, 344.451)
        bw = 70
        cl_ls = (247, 998, 8.94)
        cl_ks = (0, 0, 1.48)
        eliminations = [clearance_helper(cl_ls[i], cl_ks[i]) for i in range(len(mms))]
        partition_coefficients = partitions(logps, fus)
        self.partition_coefficients = partition_coefficients
        self.bw = bw
        blood_flows = Q_function(bw)
        volumes = V_function(bw, 3)
        self.blood_flows = blood_flows
        self.volumes = volumes
        self.eliminations = eliminations
        gravimetric_coefficients = gravimetric_model(mms)
        self.gravimetric_coefficients = gravimetric_coefficients
        self.n_of_substances = 3
        administration = placeholder_helper
        self.administration = administration
        self.odes = system_generator(volumes, blood_flows, partition_coefficients, eliminations, gravimetric_coefficients, 3, administration)


    def set_substance_and_metabolites(self, log_ps, fus, mms):
        """
        Set the substance properties that are being investigated.

        Updates partition coefficients and metabolite mass ratios based on the
        physicochemical properties and molecular masses, then regenerates the
        ODE system.

        Parameters
        ----------
        log_ps : sequence of float
            Log octanol-water partition coefficients for the substance and
            metabolites.
        fus : sequence of float
            Fractions unbound in blood, with one value for each substance.
        mms : sequence of float
            Molecular masses, with one value for each substance. The sequence
            order defines the parent substance followed by its metabolites.
        
        Notes
        -----
        This method updates the partition_coefficients and
        gravimetric_coefficients attributes and regenerates the ODE system.
        The three sequences must have the same nonzero length.
        
        Examples
        --------
        >>> m = model()
        >>> m.set_substance_and_metabolites(
        ...     (6.97, 5.33, 5.24), (0.0022448, 0.01209, 0.101),
        ...     (314.469, 330.468, 344.451))
        """

        n = len(mms)
        self.n_of_substances = n
        self.volumes = V_function(self.bw, n)
        partition_coefficients = partitions(log_ps, fus)
        self.partition_coefficients = partition_coefficients
        gravimetric_coefficients = gravimetric_model(mms)
        self.gravimetric_coefficients = gravimetric_coefficients
        self.odes = system_generator(self.volumes, self.blood_flows, partition_coefficients, self.eliminations, gravimetric_coefficients, self.n_of_substances, self.administration)


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
        Substance and metabolite properties can be set separately using
        set_substance_and_metabolites().
        
        Examples
        --------
        >>> m = model()
        >>> m.set_patient(bw=70)  # 70 kg patient
        """

        blood_flows = Q_function(bw)
        volumes = V_function(bw, self.n_of_substances)
        self.bw = bw
        self.blood_flows = blood_flows
        self.volumes = volumes
        self.odes = system_generator(volumes, blood_flows, self.partition_coefficients, self.eliminations, self.gravimetric_coefficients, self.n_of_substances, self.administration)

    def set_administration(self, roa):
        self.route_of_administration = roa
        if roa == 'inh':
            administration = inhalation_helper
        else:
            administration = placeholder_helper
        self.administration = administration
        self.odes = system_generator(self.volumes, self.blood_flows, self.partition_coefficients, self.eliminations, self.gravimetric_coefficients, self.n_of_substances, administration)
            
    def set_elimination(self, cl_ls, cl_ks):
        """
        Set the elimination parameters.

        Updates the elimination function to account for linear clearance
        from the liver and kidney. The elimination function calculates
        the outflow of substance from these tissues.

        Parameters
        ----------
        cl_ls : sequence of float
            Linear clearance values from the liver, with one value for each
            substance.
        cl_ks : sequence of float
            Linear clearance values from the kidney, with one value for each
            substance.
        
        Notes
        -----
        Each clearance function modifies compartments 5 (kidney) and 6 (liver)
        with linear kinetics: outflow = clearance * concentration. The two
        sequences must have one value per substance.
        
        Examples
        --------
        >>> m = model()
        >>> m.set_elimination((10, 10, 10), (5, 5, 5))
        """

        eliminations = [clearance_helper(cl_ls[i], cl_ks[i]) for i in range(self.n_of_substances)]
        self.eliminations = eliminations
        self.odes = system_generator(self.volumes, self.blood_flows, self.partition_coefficients, eliminations, self.gravimetric_coefficients, self.n_of_substances, self.administration)


    def simulate(self, doses, times):
        """
        Set the administration regime and simulate the system of ODEs.
        
        The `times` array must contain one more time point than the `doses`
        sequence. Each dose at index `i` is administered at `times[i]`, and
        the final entry in `times` is the last observation time. The time
        values are interpreted in minutes and must be strictly increasing.
        
        Parameters
        ----------
        doses : array-like
            A list of all doses.
        times : array-like
            A list of times at which doses are given, plus a final endpoint.
            This array must have length `len(doses) + 1` and contain
            strictly increasing values.
        route_of_administration : {'iv', 'ia', 'inh'}, optional
            Route used to administer each dose to the parent substance. 'iv'
            injects the dose into the venous blood compartment, 'ia' into the
            arterial blood compartment, and 'inh' into the lung compartment.
            Doses can only be administered to the parent substance; metabolites
            are produced by the model. The default is 'iv'.
        
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
        >>> m.set_elimination((247, 998, 8.94), (0, 0, 1.48))
        >>> t_iv, c_iv = m.simulate([0.053], [0, 24*60], route_of_administration='iv')
        >>> t_inh, c_inh = m.simulate([0.053], [0, 24*60], route_of_administration='inh')
        """
        self.t, self.c = solver(self.odes, doses, times, self.volumes, self.route_of_administration, self.n_of_substances)
        return self.t, self.c
        
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
        graph_whole_helper(self.t, self.c, name, time_unit=time_unit)

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