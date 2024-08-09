from bhpwaveformcy import (TrajectoryDataPy,
                           InspiralContainerWrapper,
                           InspiralGeneratorPy)
from bhpwave.constants import *
from bhpwave.trajectory.geodesic import (kerr_circ_geo_orbital_frequency,
                                         kerr_circ_geo_radius)
import os
import numpy as np

path_to_file = os.path.dirname(os.path.abspath(__file__))
default_path = path_to_file + "/../data/trajectory.txt"

class InspiralGenerator:
    """A class for generating quasi-circular inspirals around a rotating massive black hole.
    Once instantiated, the class can be called to generate inspiral data.

    Parameters
    ----------
    trajectory_data : TrajectoryData or None, optional
        A TrajectoryData class which holds interpolants of the relevant
        trajectory data
    """
    def __init__(self, trajectory_data = None):
        if trajectory_data is None:
            self.trajectory_data = TrajectoryDataPy(default_path, dealloc_flag=False)
        else:
            self.trajectory_data = trajectory_data.base_class
        self.inspiral_generator = InspiralGeneratorPy(self.trajectory_data)

    def check_radius(self, a, r0):
        """A utility function for checking that the orbital radius lies in the interpolation range"""
        omega_min = self.trajectory_data.min_orbital_frequency(a)
        omega_max = self.trajectory_data.max_orbital_frequency(a)
        omega = kerr_circ_geo_orbital_frequency(a, r0)
        if omega > omega_max or omega < omega_min:
            raise ValueError("Radius corresponds to a frequency that lies outside the valid domain of [{},{}]".format(omega_min, omega_max))

    def __call__(self, M, mu, a, r0, dt, T, num_threads=0):
        """Generates a quasi-circular inspiral of a point-particle around a rotating massive black hole.

        Parameters
        ----------
        M : double
            the mass of the rotating massive black hole (in solar
            masses)
        mu : double
            the mass of the smaller compact object (in solar masses)
        r0 : double
            the initial Boyer-Lindquist radius of the small body
        dt : double
            time step in seconds for sampling the trajectory
        T : double
            duration of the inspiral in years

        Returns
        -------
        Inspiral
            A class object containing inspiral data
        """
        massratio = mu/M
        dtM = dt/(M*Modot_GC1_to_S)
        TM = T*yr_MKS/(M*Modot_GC1_to_S)
        self.check_radius(a, r0)

        inspiralwrapper=self.inspiral_generator(massratio, a, r0, dtM, TM, num_threads)
        inspiral = Inspiral(inspiralwrapper)

        return inspiral

class TrajectoryData:
    """A class that holds all of the pre-computed trajectory data
    for quasi-circular inspirals around a rotating massive black hole.
    """
    def __init__(self, file_path = default_path, dealloc_flag = False):
        file_path = os.path.abspath(file_path)
        if os.path.exists(file_path):
            self.trajectory_data = TrajectoryDataPy(file_path, dealloc_flag)
        else:
            raise ValueError("File {} does not exist".format(file_path))

    def check_freq(self, a, omega):
        """A utility function for checking that the orbital frequency lies in the interpolation range"""
        omega_min = self.trajectory_data.min_orbital_frequency(a)
        omega_max = self.trajectory_data.max_orbital_frequency(a)
        if np.any(omega > omega_max) or np.any(omega < omega_min):
            raise ValueError("Radius corresponds to a frequency that lies outside the valid domain of [{},{}]".format(omega_min, omega_max))
        
    def time_to_merger(self, M, mu, a, r0):
        """Time in seconds until the system reaches the ISCO

        Parameters
        ----------
        M : double
            primary mass in solar masses
        mu : double
            secondary mass in solar masses
        a : double
            Kerr spin parameter
        r0 : double
            initial orbital radius
        """
        omega = kerr_circ_geo_orbital_frequency(a, r0)
        self.check_freq(a, omega)
        tM = self.trajectory_data.time_to_merger(a, omega)
        return tM*M*Modot_GC1_to_S/(mu/M)

    def phase_to_merger(self, M, mu, a, r0):
        """Number of orbital phase accumulated until the system reaches the ISCO

        Parameters
        ----------
        M : double
            primary mass in solar masses
        mu : double
            secondary mass in solar masses
        a : double
            Kerr spin parameter
        r0 : double
            initial orbital radius
        """
        omega = kerr_circ_geo_orbital_frequency(a, r0)
        self.check_freq(a, omega)
        phi = self.trajectory_data.phase_to_merger(a, omega)
        return phi/(mu/M)
    
    def orbital_frequency_to_merger(self, M, mu, a, T):
        """The frequency of the orbit T years before merger.

        Parameters
        ----------
        M : double
            primary mass in solar masses
        mu : double
            secondary mass in solar masses
        a : double
            Kerr spin parameter
        T : double
            years to merger
        """
        tM = -T*yr_MKS/(M*Modot_GC1_to_S/(mu/M))
        return self.trajectory_data.orbital_frequency(a, tM)
    
    def radius_to_merger(self, M, mu, a, T):
        """The radius of the orbit T years before merger.

        Parameters
        ----------
        M : double
            primary mass in solar masses
        mu : double
            secondary mass in solar masses
        a : double
            Kerr spin parameter
        T : double
            years to merger
        """
        omega = self.orbital_frequency_to_merger(M, mu, a, T)
        return kerr_circ_geo_radius(a, omega)

    def scaled_energy_flux(self, a, r0):
        """The energy flux in units :math:`G = c = 1` and scaled by the mass ratio

        Parameters
        ----------
        a : double
            Kerr spin parameter
        r0 : double
            initial orbital radius
        """
        MIN_ARRAY_LENGTH = 10
        omega = kerr_circ_geo_orbital_frequency(a, r0)
        # self.check_freq(a, omega)
        if isinstance(a, np.ndarray) and isinstance(omega, np.ndarray):
            assert a.shape == omega.shape
            dim = a.shape
            if len(dim) == 1:
                if a.shape[0] > MIN_ARRAY_LENGTH:
                    return self.trajectory_data.flux_parallel(a, omega)
                else:
                    return np.array([self.trajectory_data.flux(a[i], omega[i]) for i in range(a.shape[0])])
            elif len(dim) == 2:
                arr = self.trajectory_data.flux_parallel(a.flatten(), omega.flatten())
                return arr.reshape(dim)
            else:
                return ValueError("Can only handle numpy arrays in up to 2-dimensions.")
        elif isinstance(omega, np.ndarray):
            if omega.shape[0] > MIN_ARRAY_LENGTH:
                return self.trajectory_data.flux_parallel_frequency(a, omega)
            else:
                return np.array([self.trajectory_data.flux(a[i], omega[i]) for i in range(omega.shape[0])])
        elif isinstance(a, np.ndarray):
            if a.shape[0] > MIN_ARRAY_LENGTH:
                return self.trajectory_data.flux_parallel_spin(a, omega)
            else:
                return np.array([self.trajectory_data.flux(a[i], omega[i]) for i in range(a.shape[0])])
        else:
            return self.trajectory_data.flux(a, omega)

    @property
    def base_class(self):
        """Returns the base Cython class"""
        return self.trajectory_data

class Inspiral:
    """A class that holds inspiral output from the InspiralGenerator."""
    def __init__(self, inspiral_wrapper):
        self.inspiral_data = inspiral_wrapper

    @property
    def size(self):
        """Size of the arrays within the Inspiral class. Equivalent to the number of time steps."""
        return self.inspiral_data.size
    
    @property
    def time(self):
        """Evolution of time"""
        return self.inspiral_data.time

    @property
    def frequency(self):
        """Evolution of the orbital frequency"""
        return self.inspiral_data.frequency

    @property
    def radius(self):
        """Evolution of the orbital radius"""
        return self.inspiral_data.radius
    
    @property
    def phase(self):
        """Evolution of the orbital phase"""
        return self.inspiral_data.phase

    @property
    def spin(self):
        """Black hole spin of the system"""
        return self.inspiral_data.spin
    
    @property
    def a(self):
        """Black hole spin of the system"""
        return self.inspiral_data.a
    
    @property
    def massratio(self):
        """Mass ratio of the system"""
        return self.inspiral_data.massratio
    
    @property
    def initialradius(self):
        """Initial orbital radius of the inspiral"""
        return self.inspiral_data.initialradius

    @property
    def initialfrequency(self):
        """Initial orbital frequency of the inspiral"""
        return self.inspiral_data.initialfrequency

    @property
    def iscofrequency(self):
        """Frequency of the innermost stable circular orbit for this spacetime"""
        return self.inspiral_data.iscofrequency

    @property
    def dt(self):
        """Size of the time steps in the inspiral"""
        return self.inspiral_data.dt

    @property
    def base_class(self):
        """Returns the base Cython class"""
        return self.inspiral_data