from bhpwaveformcy import (TrajectoryDataPy,
                           InspiralContainerWrapper,
                           InspiralGeneratorPy)
from bhpwave.constants import *
from bhpwave.trajectory.geodesic import kerr_circ_geo_orbital_frequency, kerr_isco_frequency, kerr_circ_geo_radius
import os

from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
import numpy as np

path_to_file = os.path.dirname(os.path.abspath(__file__))
default_path = path_to_file + "/../data/trajectory.txt"

class InspiralGenerator:
    """
    A class for generating quasi-circular inspirals around a rotating massive black hole.
    Once instantiated, the class can be called to generate inspiral data.

    :param trajectory_data: A TrajectoryData class which holds interpolants of the relevant trajectory data
    :type trajectory_data: TrajectoryData or None, optional
    """
    def __init__(self, trajectory_data = None):
        if trajectory_data is None:
            self.trajectory_data = TrajectoryDataPy(default_path, dealloc_flag=False)
        else:
            self.trajectory_data = trajectory_data.base_class
        self.inspiral_generator = InspiralGeneratorPy(self.trajectory_data)

    def check_radius(self, a, r0):
        """
        A utility function for checking that the orbital radius lies in the interpolation range
        """
        omega_min = self.trajectory_data.min_orbital_frequency(a)
        omega_max = self.trajectory_data.max_orbital_frequency(a)
        omega = kerr_circ_geo_orbital_frequency(a, r0)
        if omega > omega_max or omega < omega_min:
            raise ValueError("Radius corresponds to a frequency that lies outside the valid domain of [{},{}]".format(omega_min, omega_max))

    def __call__(self, M, mu, a, r0, dt, T, num_threads=0):
        """
        Generates a quasi-circular inspiral of a point-particle around a rotating massive black hole.

        :param M: the mass of the rotating massive black hole (in solar masses)
        :type M: double
        :param mu: the mass of the smaller compact object (in solar masses)
        :type mu: double
        :param r0: the initial Boyer-Lindquist radius of the small body
        :type r0: double
        :param dt: time step in seconds for sampling the trajectory
        :type dt: double
        :param T: duration of the inspiral in years
        :type T: double

        :return: A class object containing inspiral data
        :rtype: Inspiral
        """
        massratio = mu/M
        dtM = dt/(M*Modot_GC1_to_S)
        TM = T*yr_MKS/(M*Modot_GC1_to_S)
        self.check_radius(a, r0)

        inspiralwrapper=self.inspiral_generator(massratio, a, r0, dtM, TM, num_threads)
        inspiral = Inspiral(inspiralwrapper)

        return inspiral
    
class IntegrateInspiralGeneratorBase:
    """
    A class for generating quasi-circular inspirals around a rotating massive black hole.
    Once instantiated, the class can be called to generate inspiral data.

    """
    def __init__(self, Edot_norm = None):
        self.Edot_norm = Edot_norm

    def isco_radius(self, a):
        return a

    def omega_of_a_alpha(self, a, alpha):
        return 0.

    def chi_of_spin(self, a):
        return (1. - a)**(1./3.)
    
    def alpha_of_omega(self, omega, oISCO):
        return 0.
    
    def omega_of_a_r(self, a, r):
        return 0.
    
    def omega_alpha_derivative(self, omega, oISCO):
        return 0.
    
    def energy_omega_derivative(self, a, omega):
        return 0.
    
    def pn_flux_noprefactor(self, omega):
        return 0.

    def _dJdTIntegrate(self, t, alphaVec, a, oISCO, *args):
        alpha = alphaVec[0]
        chi = self.chi_of_spin(a)
        omega = self.omega_of_a_alpha(a, alpha)
        dOmega_dAlpha = self.omega_alpha_derivative(omega, oISCO)
        dE_dOmega = self.energy_omega_derivative(a, omega)
        Edot = self.Edot_func(chi, alpha)*self.pn_flux_noprefactor(omega)
        return np.array([1./(dE_dOmega*dOmega_dAlpha/Edot), omega])

    def integrate_eom(self, M, mu, a, r0, dt, T, *args):
        omega0 = self.omega_of_a_r(a, r0)
        rISCO = self.isco_radius(a)
        oISCO = self.omega_of_a_r(a, rISCO)
        alpha0 = self.alpha_of_omega(omega0, oISCO)

        t_array = np.linspace(0., T, int(T/dt + 1))
        insp = solve_ivp(self._dJdTIntegrate, [t_array[0], t_array[-1]], [alpha0, 0.], args=(a, oISCO, *args), method='DOP853', t_eval=t_array, rtol=1.e-13, atol = 1.e-12)
        alpha = insp.y[0]
        Phi = insp.y[1]
        inspiral = InspiralContainerWrapper(alpha.size)
        inspiral.set_initial_conditions(a, mu/M, r0, dt)
        for i in range(inspiral.timesteps):
            inspiral.set_time_step(i, alpha[i], Phi[i])

        return inspiral

    def __call__(self, M, mu, a, r0, dt, T, *args):
        """
        Generates a quasi-circular inspiral of a point-particle around a rotating massive black hole.

        :param M: the mass of the rotating massive black hole (in solar masses)
        :type M: double
        :param mu: the mass of the smaller compact object (in solar masses)
        :type mu: double
        :param r0: the initial Boyer-Lindquist radius of the small body
        :type r0: double
        :param dt: time step in seconds for sampling the trajectory
        :type dt: double
        :param T: duration of the inspiral in years
        :type T: double

        :return: A class object containing inspiral data
        :rtype: Inspiral
        """
        massratio = mu/M
        dtM = dt/(M*Modot_GC1_to_S)
        TM = T*yr_MKS/(M*Modot_GC1_to_S)

        inspiralwrapper = self.integrate_eom(M, mu, a, r0, dtM, TM, *args)
        inspiral = Inspiral(inspiralwrapper)

        return inspiral

class TrajectoryData:
    """
    A class that holds all of the pre-computed trajectory data
    for quasi-circular inspirals around a rotating massive black hole.
    """
    def __init__(self, file_path = default_path, dealloc_flag = False):
        file_path = os.path.abspath(file_path)
        if os.path.exists(file_path):
            self.trajectory_data = TrajectoryDataPy(file_path, dealloc_flag)
        else:
            raise ValueError("File {} does not exist".format(file_path))

    def check_freq(self, a, omega):
        """
        A utility function for checking that the orbital frequency lies in the interpolation range
        """
        omega_min = self.trajectory_data.min_orbital_frequency(a)
        omega_max = self.trajectory_data.max_orbital_frequency(a)
        if omega > omega_max or omega < omega_min:
            raise ValueError("Radius corresponds to a frequency that lies outside the valid domain of [{},{}]".format(omega_min, omega_max))
        
    def time_to_merger(self, M, mu, a, r0):
        """
        Time in seconds until the system reaches the ISCO

        :param M: primary mass in solar masses
        :type M: double
        :param mu: secondary mass in solar masses
        :type mu: double
        :param a: Kerr spin parameter
        :type a: double
        :param r0: initial orbital radius
        :type r0: double
        """
        omega = kerr_circ_geo_orbital_frequency(a, r0)
        self.check_freq(a, omega)
        tM = self.trajectory_data.time_to_merger(a, omega)
        return tM*M*Modot_GC1_to_S/(mu/M)

    def phase_to_merger(self, M, mu, a, r0):
        """
        Number of orbital phase accumulated until the system reaches the ISCO

        :param M: primary mass in solar masses
        :type M: double
        :param mu: secondary mass in solar masses
        :type mu: double
        :param a: Kerr spin parameter
        :type a: double
        :param r0: initial orbital radius
        :type r0: double
        """
        omega = kerr_circ_geo_orbital_frequency(a, r0)
        self.check_freq(a, omega)
        phi = self.trajectory_data.phase_to_merger(a, omega)
        return phi/(mu/M)

    def scaled_energy_flux(self, a, r0):
        """
        The energy flux in units G = c = 1 and scaled by the mass ratio

        :param a: Kerr spin parameter
        :type a: double
        :param r0: initial orbital radius
        :type r0: double
        """
        omega = kerr_circ_geo_orbital_frequency(a, r0)
        self.check_freq(a, omega)
        return self.trajectory_data.flux(a, omega)

    @property
    def base_class(self):
        """
        Returns the base Cython class
        """
        return self.trajectory_data

class Inspiral:
    """
    A class that holds inspiral output from the InspiralGenerator.
    """
    def __init__(self, inspiral_wrapper):
        self.data = inspiral_wrapper

    @property
    def size(self):
        """
        Size of the arrays within the Inspiral class. Equivalent to the number of time steps.
        """
        return self.data.size
    
    @property
    def time(self):
        """
        Evolution of time
        """
        return self.data.time

    @property
    def frequency(self):
        """
        Evolution of the orbital frequency
        """
        return self.data.frequency

    @property
    def radius(self):
        """
        Evolution of the orbital radius
        """
        return self.data.radius
    
    @property
    def phase(self):
        """
        Evolution of the orbital phase
        """
        return self.data.phase

    @property
    def spin(self):
        """
        Black hole spin of the system
        """
        return self.data.spin
    
    @property
    def a(self):
        """
        Black hole spin of the system
        """
        return self.data.a
    
    @property
    def massratio(self):
        """
        Mass ratio of the system
        """
        return self.data.massratio
    
    @property
    def initialradius(self):
        """
        Initial orbital radius of the inspiral
        """
        return self.data.initialradius

    @property
    def initialfrequency(self):
        """
        Initial orbital frequency of the inspiral
        """
        return self.data.initialfrequency

    @property
    def iscofrequency(self):
        """
        Frequency of the innermost stable circular orbit for this spacetime
        """
        return self.data.iscofrequency

    @property
    def dt(self):
        """
        Size of the time steps in the inspiral
        """
        return self.data.dt

    @property
    def base_class(self):
        """
        Returns the base Cython class
        """
        return self.data
    
    @property
    def base(self):
        """
        Returns the base Cython class
        """
        return self.data