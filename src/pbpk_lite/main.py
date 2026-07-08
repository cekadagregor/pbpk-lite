"""
Main module for PBPK (Physiologically Based Pharmacokinetic) modeling.

This module provides the core model class for setting up and solving
PBPK models with customizable substance properties, patient parameters,
and elimination pathways.
"""

from .distribution_params import V_function, Q_function, partition_model
from .graph import graph_whole_helper, graph_venous_helper, graph_compartments_helper
from .solve import solver
from .ode_system import system_generator
from .elimination import clearance_helper, michaelis_menten_helper

class model:
    """
    Physiologically Based Pharmacokinetic (PBPK) model class.
    
    This class manages the setup and solution of PBPK models, handling
    substance properties, patient parameters, and drug elimination pathways.
    It integrates partition coefficients, blood flows, tissue volumes, and
    elimination kinetics to solve the system of ordinary differential equations.
    The implementation assumes a consistent unit system based on mass units
    for doses and amounts, mL for volumes, and minutes for time.
    
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
    >>> doses = [0.053]
    >>> times = [0, 24*60]
    >>> t, c = m.simulate(doses, times)
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
        Elimination is initialized with a placeholder clearance value for 
        cl_l and a zero value for cl_k which can be updated using 
        set_elimination().
        """
        # placeholder parameters
        logp = 6.97
        fu = 0.0022448
        bw = 70
        cl_l = 748.6643482986731
        cl_k = 0
        self.partition_coefficients = partition_model(logp, fu)
        self.blood_flows = Q_function(bw)
        self.volumes = V_function(bw)
        self.elimination = clearance_helper(cl_l, cl_k)
        self.odes = system_generator(self.volumes, self.blood_flows, self.partition_coefficients, self.elimination)

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
        The elimination function modifies compartments 5 (kidney) and 6 (liver)
        with linear kinetics: outflow = clearance * concentration.
        After setting elimination parameters, the ODE system is regenerated.
        
        Examples
        --------
        >>> m = model()
        >>> m.set_elimination(cl_l=10, cl_k=5)
        """

        elimination = clearance_helper(cl_l, cl_k)
        self.odes = system_generator(self.volumes, self.blood_flows, self.partition_coefficients, elimination)
    
    def set_michaelis_menten_elimination(self, vmax, km):
        """
        Set the elimination parameters using Michaelis-Menten kinetics.

        Updates the elimination function to account for nonlinear clearance
        from the liver based on Michaelis-Menten kinetics. The elimination
        function calculates the outflow of substance from the liver compartment.

        Parameters
        ----------
        vmax : float
            Maximum rate of elimination (Vmax) from the liver.
        km : float
            Michaelis constant (Km) representing the concentration at which
            the elimination rate is half of Vmax.
        
        Notes
        -----
        The elimination function modifies compartment 6 (liver) with nonlinear
        kinetics: outflow = Vmax * concentration / (Km + concentration).
        After setting Michaelis-Menten elimination parameters, the ODE system
        is regenerated.
        
        Examples
        --------
        >>> m = model()
        >>> m.set_michaelis_menten_elimination(vmax=10, km=5)
        """
        
        elimination = michaelis_menten_helper(vmax, km)
        self.odes = system_generator(self.volumes, self.blood_flows, self.partition_coefficients, elimination)

    def simulate(self, doses, times, route_of_administration='iv'):
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
            Route used to administer each dose. 'iv' injects the dose into
            the venous blood compartment, 'ia' into the arterial blood
            compartment, and 'inh' into the lung compartment. The default is
            'iv'.
        
        Returns
        -------
        t : ndarray
            Times at which concentrations are calculated.
        c : ndarray
            Concentrations in each tissue at the given times.
        
        Examples
        --------
        >>> m = model()
        >>> m.set_substance(log_p=6.97, fu=0.0022448)
        >>> m.set_patient(bw=70)
        >>> t_iv, c_iv = m.simulate([0.053], [0, 24*60], route_of_administration='iv')
        >>> t_inh, c_inh = m.simulate([0.053], [0, 24*60], route_of_administration='inh')
        """

        self.t, self.c = solver(self.odes, doses, times, self.volumes, route_of_administration)
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