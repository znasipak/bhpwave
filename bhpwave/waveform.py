import numpy as np
from bhpwaveformcy import (WaveformHarmonicGeneratorPy,
                           WaveformFourierGeneratorPy,
                           WaveformGeneratorPy,
                           TrajectoryDataPy,
                           InspiralGeneratorPy,
                           HarmonicAmplitudesPy)

from bhpwave.trajectory.geodesic import A_MAX
from bhpwave.trajectory.inspiral import TrajectoryData
from bhpwave.harmonics.amplitudes import HarmonicAmplitudes

import os
from bhpwave.constants import *
import warnings

path_to_file = os.path.dirname(os.path.abspath(__file__))
traj_path = path_to_file + "/data/trajectory.txt"
amplitude_path = path_to_file + "/data/circ_data"

try:
    CPU_MAX = len(os.sched_getaffinity(0))
except:
    CPU_MAX = os.cpu_count()

def check_source_frame_parameters(a, x0, theta, phi, Phi_phi0):
    assert np.abs(x0) == 1
    assert np.abs(a) <= A_MAX

    if x0 < 0:
        a = -a
        x0 = -x0
        theta = np.pi - theta
        phi = -phi
        Phi_phi0 = -Phi_phi0

    return a, x0, theta, phi, Phi_phi0

def check_ssb_frame_parameters(a, x0, qK, phiK, Phi_phi0):
    assert np.abs(x0) == 1
    assert np.abs(a) <= A_MAX

    if x0 < 0:
        a = -a
        x0 = -x0
        qK = np.pi - qK
        phiK = np.pi + phiK
        Phi_phi0 = np.pi + Phi_phi0

    return a, x0, qK, phiK, Phi_phi0

class KerrCircularWaveformCustomTrajectory:
    def __init__(self, trajectory_class, harmonic_data, waveform_kwargs):
        self.waveform_generator = WaveformHarmonicGeneratorPy(harmonic_data, waveform_kwargs)
        self.inspiral_generator = trajectory_class

    def select_modes(self, M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs):
        inspiral = self.inspiral_generator(M, mu, a, r0, dt, T, *args, **kwargs)
        return self.waveform_generator.select_modes(inspiral.data, qS, phiS, qK, phiK, **kwargs)
    
    def waveform_harmonics(self, l, m, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs):
        inspiral = self.inspiral_generator(M, mu, a, r0, dt, T, *args, **kwargs)
        return self.waveform_generator.waveform_harmonics(l, m, mu, inspiral.data, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, **kwargs)
    
    def waveform(self, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs):
        inspiral = self.inspiral_generator(M, mu, a, r0, dt, T, *args, **kwargs)
        return self.waveform_generator.waveform(mu, inspiral.data, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, **kwargs)
    
    def waveform_harmonics_grid(self, l, m, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs):
        inspiral = self.inspiral_generator(M, mu, a, r0, dt, T, *args, **kwargs)
        return self.waveform_generator.waveform_harmonics_grid(l, m, mu, inspiral.data, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, **kwargs)
    
    def waveform_select_harmonics_grid(self, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs):
        inspiral = self.inspiral_generator(M, mu, a, r0, dt, T, *args, **kwargs)
        return self.waveform_generator.waveform_select_harmonics_grid(mu, inspiral.data, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, **kwargs)
    
    def waveform_harmonics_phase_amplitude(self, l, m, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs):
        inspiral = self.inspiral_generator(M, mu, a, r0, dt, T, *args, **kwargs)
        return self.waveform_generator.waveform_harmonics_phase_amplitude(l, m, mu, inspiral.data, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, **kwargs)
    
    def waveform_select_harmonics_phase_amplitude(self, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs):
        inspiral = self.inspiral_generator(M, mu, a, r0, dt, T, *args, **kwargs)
        return self.waveform_generator.waveform_select_harmonics_phase_amplitude(mu, inspiral.data, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, **kwargs)
    
    def waveform_source_frame_harmonics(self, l, m, M, mu, a, r0, theta, phi, Phi_phi0, dt, T, *args, **kwargs):
        inspiral = self.inspiral_generator(M, mu, a, r0, dt, T, *args, **kwargs)
        return self.waveform_generator.waveform_harmonics_phase_amplitude(l, m, inspiral.data, theta, phi, Phi_phi0, dt, T, **kwargs)
    
    def waveform_source_frame(self, M, mu, a, r0, theta, phi, Phi_phi0, dt, T, *args, **kwargs):
        inspiral = self.inspiral_generator(M, mu, a, r0, dt, T, *args, **kwargs)
        return self.waveform_generator.waveform_select_harmonics_phase_amplitude(inspiral.data, theta, phi, Phi_phi0, dt, T, **kwargs)

class KerrCircularWaveform:
    """
    Class that generates the gravitational waveform produced by an extreme-mass-ratio inspiral
    using the adiabatic approximation from black hole perturbation theory and the self-force formalism.
    By default, the waveform is generated in the solar system barycenter frame.

    Waveform generation is limited to quasi-circular inspirals in Kerr spacetime, but the generator mirrors the generic
    parametrization used in other EMRI waveform generators (e.g., https://bhptoolkit.org/FastEMRIWaveforms/html/user/main.html)

    :param trajectory_data: a TrajectoryData class which holds interpolants of the relevant trajectory data
    :type trajectory_data: TrajectoryData or None, optional
    :param harmonic_data: a HarmonicAmplitudes class which holds interpolants of the harmonic mode amplitudes
    :type harmonic_data: HarmonicAmplitudes or None, optional
    :param num_threads: the number of threads used to evaluate the waveform
    :type num_threads: int or None, optional
    """
    def __init__(self, trajectory_data=None, harmonic_data=None, num_threads=None):
        self.custom_inspiral_flag = False
        if num_threads is None:
            num_threads = CPU_MAX
        if trajectory_data is None:
            self.trajectory_data = TrajectoryDataPy(filename=traj_path, dealloc_flag=False)
        else:
            if isinstance(trajectory_data, TrajectoryData):
                self.trajectory_data = trajectory_data.base_class
            else:
                self.custom_inspiral_flag = True
                self.trajectory_data = trajectory_data
            
        if harmonic_data is None:
            self.harmonic_data = HarmonicAmplitudesPy(filebase=amplitude_path, dealloc_flag=False)
        else:
            if isinstance(harmonic_data, HarmonicAmplitudes):
                self.harmonic_data = harmonic_data.base_class
            else:
                raise ValueError("Wrong type for harmonic_data. Must be an instance of the HarmonicAmplitudes class.")

        waveform_kwargs = {
            "num_threads": num_threads
        }

        if self.custom_inspiral_flag:
            self.waveform_generator = KerrCircularWaveformCustomTrajectory(self.trajectory_data, self.harmonic_data, waveform_kwargs=waveform_kwargs)
        else:
            self.waveform_generator = WaveformGeneratorPy(self.trajectory_data, self.harmonic_data, waveform_kwargs=waveform_kwargs)

    def select_modes(self, M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs):
        """
        Selects the harmonic modes that are used for calculating the waveform

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param r0: initial orbital separation of the two objects
        :type r0: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param T: Duration of the waveform in years
        :type T: double, optional

        :rtype: 1d-array[tuples(doubles)]
        """

        return self.waveform_generator.select_modes(M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs)

    def __call__(self, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs):
        """
        Calculate the complex gravitational wave strain

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param r0: initial orbital separation of the two objects
        :type r0: double
        :param dist: luminosity distance to the source in Gpc
        :type dist: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param dt: Spacing of time samples in seconds
        :type dt: double, optional
        :param T: Duration of the waveform in years
        :type T: double, optional
        
        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        :param include_negative_m: True returns the sum of the positive and negative m-modes for each mode in select_modes
        :type include_negative_m: bool, optional
        
        :rtype: 1d-array[complex] or list[two 1d-arrays[double]]

        """

        if "T" in kwargs.keys():
            T = kwargs["T"]
        if "dt" in kwargs.keys():
            dt = kwargs["dt"]
        
        include_negative_m = True
        if "include_negative_m" in kwargs.keys():
            include_negative_m = kwargs["include_negative_m"]

        if "select_modes" in kwargs.keys():
            lmodes = []
            mmodes = []
            if include_negative_m:
                for mode in kwargs["select_modes"]:
                    if mode[1] > 0: # if include_negative_m is True then only keep positive m
                        lmodes.append(mode[0])
                        mmodes.append(mode[1])
                    else:
                        warnings.warn("Warning: Only keeping modes in select_modes with m > 0. Set include_negative_m = False to keep m < 0 modes.")
            else: # if include_negative_m is False then keep all m
                for mode in kwargs["select_modes"]:
                    lmodes.append(mode[0])
                    mmodes.append(mode[1])
            l = np.ascontiguousarray(lmodes, dtype=np.intc)
            m = np.ascontiguousarray(mmodes, dtype=np.intc)
            h = self.waveform_generator.waveform_harmonics(l, m, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs)
        else:
            h = self.waveform_generator.waveform(M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs)
        return h
    
    def harmonics(self, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs):
        """
        Calculate the spin-weighted spherical harmonic modes of the gravitational wave strain

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param r0: initial orbital separation of the two objects
        :type r0: double
        :param dist: luminosity distance to the source in Gpc
        :type dist: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param dt: Spacing of time samples in seconds
        :type dt: double, optional
        :param T: Duration of the waveform in years
        :type T: double, optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        :param include_negative_m: True returns the sum of the positive and negative m-modes for each mode in select_modes
        :type include_negative_m: bool, optional
        
        :rtype: 2d-array[complex] or list[two 2d-arrays[double]]

        """

        if "T" in kwargs.keys():
            T = kwargs["T"]
        if "dt" in kwargs.keys():
            dt = kwargs["dt"]

        include_negative_m = True
        if "include_negative_m" in kwargs.keys():
            include_negative_m = kwargs["include_negative_m"]

        if "select_modes" in kwargs.keys():
            lmodes = []
            mmodes = []
            if include_negative_m:
                for mode in kwargs["select_modes"]:
                    if mode[1] > 0: # if include_negative_m is True then only keep positive m
                        lmodes.append(mode[0])
                        mmodes.append(mode[1])
                    else:
                        warnings.warn("Warning: Only keeping modes in select_modes with m > 0. Set include_negative_m = False to keep m < 0 modes.")
            else: # if include_negative_m is False then keep all m
                for mode in kwargs["select_modes"]:
                    lmodes.append(mode[0])
                    mmodes.append(mode[1])
            l = np.ascontiguousarray(lmodes, dtype=np.intc)
            m = np.ascontiguousarray(mmodes, dtype=np.intc)
            h = self.waveform_generator.waveform_harmonics_grid(l, m, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs)
        else:
            h = self.waveform_generator.waveform_select_harmonics_grid(M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs)
        return h
    
    def harmonics_data(self, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs):
        """
        Calculate the amplitudes and phases of the spin-weighted spherical harmonic (l, m)-modes of the gravitational wave strain. The
        function returns mode data for both positive and negative m-modes. Given a list of modes [(l_1, m_1), (l_2, m_2), ... , (l_N, m_N)]
        the function returns data in the order [(l_1, m_1), (l_1, -m_1), (l_2, m_2), (l_2, -m_2), ..., (l_N, m_N), (l_N, -m_N)].

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param r0: initial orbital separation of the two objects
        :type r0: double
        :param dist: luminosity distance to the source in Gpc
        :type dist: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param dt: Spacing of time samples in seconds
        :type dt: double, optional
        :param T: Duration of the waveform in years
        :type T: double, optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation. Only positive m-modes are used.
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        
        :rtype: 2d-array[complex] or list[two 2d-arrays[double]]

        """

        if "T" in kwargs.keys():
            T = kwargs["T"]
        if "dt" in kwargs.keys():
            dt = kwargs["dt"]

        if "select_modes" in kwargs.keys():
            lmodes = []
            mmodes = []
            for mode in kwargs["select_modes"]:
                if mode[1] > 0: # only keep positive m
                    lmodes.append(mode[0])
                    mmodes.append(mode[1])
                else:
                    warnings.warn("Warning: Only keeping modes in select_modes with m > 0.")
            l = np.ascontiguousarray(lmodes, dtype=np.intc)
            m = np.ascontiguousarray(mmodes, dtype=np.intc)
            h = self.waveform_generator.waveform_harmonics_phase_amplitude(l, m, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs)
        else:
            h = self.waveform_generator.waveform_select_harmonics_phase_amplitude(M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs)
        return h
    
    def source_frame(self, M, mu, a, r0, theta, phi, Phi_phi0, dt=10., T=1., *args, **kwargs):
        """
        Calculate the scaled gravitational wave strain :math:`r\\times h/\\mu` in the source frame

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param r0: initial radial separation
        :type r0: double
        :param theta: polar angle of the observor with respect to the Kerr spin vector
        :type theta: double
        :param phi: azimuthal angle of the observor with respect to the Kerr spin vector
        :type phi: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param dt: Spacing of time samples in seconds. Default is 10 seconds.
        :type dt: double, optional
        :param T: Duration of the waveform in years. Default is 1 year.
        :type T: double, optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        :param include_negative_m: True returns the sum of the positive and negative m-modes for each mode in select_modes
        :type include_negative_m: bool, optional
        
        :rtype: 1d-array[complex] or list[two 1d-arrays[double]]

        """

        if "T" in kwargs.keys():
            T = kwargs["T"]
        if "dt" in kwargs.keys():
            dt = kwargs["dt"]

        include_negative_m = True
        if "include_negative_m" in kwargs.keys():
            include_negative_m = kwargs["include_negative_m"]

        if "select_modes" in kwargs.keys():
            lmodes = []
            mmodes = []
            if include_negative_m:
                for mode in kwargs["select_modes"]:
                    if mode[1] > 0: # if include_negative_m is True then only keep positive m
                        lmodes.append(mode[0])
                        mmodes.append(mode[1])
                    else:
                        warnings.warn("Warning: Only keeping modes in select_modes with m > 0. Set include_negative_m = False to keep m < 0 modes.")
            else: # if include_negative_m is False then keep all m
                for mode in kwargs["select_modes"]:
                    lmodes.append(mode[0])
                    mmodes.append(mode[1])
            l = np.ascontiguousarray(lmodes, dtype=np.intc)
            m = np.ascontiguousarray(mmodes, dtype=np.intc)
            h = self.waveform_generator.waveform_harmonics_source_frame(l, m, M, mu, a, r0, theta, phi, Phi_phi0, dt, T, *args, **kwargs)
        else:
            h = self.waveform_generator.waveform_source_frame(M, mu, a, r0, theta, phi, Phi_phi0, dt, T, *args, **kwargs)
        return h

class KerrWaveform(KerrCircularWaveform):
    """
    Class that generates the gravitational waveform produced by an extreme-mass-ratio inspiral
    using the adiabatic approximation from black hole perturbation theory and the self-force formalism.
    By default, the waveform is generated in the solar system barycenter frame.
    
    Waveform generation is limited to quasi-circular inspirals in Kerr spacetime, but the generator mirrors the generic
    parametrization used in other EMRI waveform generators (e.g., https://bhptoolkit.org/FastEMRIWaveforms/html/user/main.html)

    :param trajectory_data: a TrajectoryData class which holds interpolants of the relevant trajectory data
    :type trajectory_data: TrajectoryData or None, optional
    :param harmonic_data: a HarmonicAmplitudes class which holds interpolants of the harmonic mode amplitudes
    :type harmonic_data: HarmonicAmplitudes or None, optional
    :param num_threads: the number of threads used to evaluate the waveform
    :type num_threads: int or None, optional  
    """
    def __call__(self, M, mu, a, p0, e0, x0, dist, qS, phiS, qK, phiK, Phi_phi0, Phi_r0, Phi_theta0, dt=10., T=1., *args, **kwargs):
        """
        Calculate the gravitational wave strain

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param p0: initial semi-latus rectum
        :type p0: double
        :param e0: initial orbital eccentricity
        :type e0: double
        :param x0: intial cosine of the orbital inclination
        :type x0: double
        :param dist: luminosity distance to the source in Gpc
        :type dist: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param Phi_r0: Phase describing the initial radial position and velocity of the small compact object
        :type Phi_r0: double
        :param Phi_theta0: Phase describing the initial polar position and velocity of the small compact object
        :type Phi_theta0: double
        :param dt: Spacing of time samples in seconds
        :type dt: double, optional
        :param T: Duration of the waveform in years
        :type T: double, optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        :param include_negative_m: True returns the sum of the positive and negative m-modes for each mode in select_modes
        :type include_negative_m: bool, optional
        
        :rtype: 1d-array[complex] or list[two 1d-arrays[double]]

        """
        a, x0, qK, phiK, Phi_phi0 = check_ssb_frame_parameters(a, x0, qK, phiK, Phi_phi0)
        return super().__call__(M, mu, a, p0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs)
    
    def harmonics(self, M, mu, a, p0, e0, x0, dist, qS, phiS, qK, phiK, Phi_phi0, Phi_r0, Phi_theta0, dt=10., T=1., *args, **kwargs):
        """
        Calculate the spin-weighted spherical harmonic modes of the gravitational wave strain

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param p0: initial semi-latus rectum
        :type p0: double
        :param e0: initial orbital eccentricity
        :type e0: double
        :param x0: intial cosine of the orbital inclination
        :type x0: double
        :param dist: luminosity distance to the source in Gpc
        :type dist: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param Phi_r0: Phase describing the initial radial position and velocity of the small compact object
        :type Phi_r0: double
        :param Phi_theta0: Phase describing the initial polar position and velocity of the small compact object
        :type Phi_theta0: double
        :param dt: Spacing of time samples in seconds
        :type dt: double, optional
        :param T: Duration of the waveform in years
        :type T: double, optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        :param include_negative_m: True returns the sum of the positive and negative m-modes for each mode in select_modes
        :type include_negative_m: bool, optional
        
        :rtype: 2d-array[complex] or list[two 2d-arrays[double]]

        """
        a, x0, qK, phiK, Phi_phi0 = check_ssb_frame_parameters(a, x0, qK, phiK, Phi_phi0)
        return super().harmonics(M, mu, a, p0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs)
    
    def harmonics_data(self, M, mu, a, p0, e0, x0, dist, qS, phiS, qK, phiK, Phi_phi0, Phi_r0, Phi_theta0, dt=10., T=1., *args, **kwargs):
        """
        Calculate the amplitudes and phases of the spin-weighted spherical harmonic (l, m)-modes of the gravitational wave strain. The
        function returns mode data for both positive and negative m-modes. Given a list of modes [(l_1, m_1), (l_2, m_2), ... , (l_N, m_N)]
        the function returns data in the order [(l_1, m_1), (l_1, -m_1), (l_2, m_2), (l_2, -m_2), ..., (l_N, m_N), (l_N, -m_N)].

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param p0: initial semi-latus rectum
        :type p0: double
        :param e0: initial orbital eccentricity
        :type e0: double
        :param x0: intial cosine of the orbital inclination
        :type x0: double
        :param dist: luminosity distance to the source in Gpc
        :type dist: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param Phi_r0: Phase describing the initial radial position and velocity of the small compact object
        :type Phi_r0: double
        :param Phi_theta0: Phase describing the initial polar position and velocity of the small compact object
        :type Phi_theta0: double
        :param dt: Spacing of time samples in seconds
        :type dt: double, optional
        :param T: Duration of the waveform in years
        :type T: double, optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation. Only positive m-modes are used.
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        
        :rtype: 2d-array[complex] or list[two 2d-arrays[double]]

        """
        a, x0, qK, phiK, Phi_phi0 = check_ssb_frame_parameters(a, x0, qK, phiK, Phi_phi0)
        return super().harmonics_data(M, mu, a, p0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, *args, **kwargs)

    
class KerrCircularFrequencyWaveform:
    """
    Class that generates the gravitational waveform produced by an extreme-mass-ratio inspiral
    in the frequency domain using the adiabatic approximation from black hole perturbation theory and the self-force formalism.
    By default, the waveform is generated in the solar system barycenter frame.

    Waveform generation is limited to quasi-circular inspirals in Kerr spacetime, but the generator mirrors the generic
    parametrization used in other EMRI waveform generators (e.g., https://bhptoolkit.org/FastEMRIWaveforms/html/user/main.html)

    :param trajectory_data: a TrajectoryData class which holds interpolants of the relevant trajectory data
    :type trajectory_data: TrajectoryData or None, optional
    :param harmonic_data: a HarmonicAmplitudes class which holds interpolants of the harmonic mode amplitudes
    :type harmonic_data: HarmonicAmplitudes or None, optional
    :param num_threads: the number of threads used to evaluate the waveform
    :type num_threads: int or None, optional
    """
    def __init__(self, trajectory_data=None, harmonic_data=None, num_threads=None):
        if num_threads is None:
            num_threads = CPU_MAX
        if trajectory_data is None:
            self.trajectory_data = TrajectoryDataPy(filename=traj_path, dealloc_flag=False)
        else:
            self.trajectory_data = trajectory_data.base_class
        if harmonic_data is None:
            self.harmonic_data = HarmonicAmplitudesPy(filebase=amplitude_path, dealloc_flag=False)
        else:
            self.harmonic_data = harmonic_data.base_class

        waveform_kwargs = {
            "num_threads": num_threads
        }
        self.waveform_generator = WaveformFourierGeneratorPy(self.trajectory_data, self.harmonic_data, waveform_kwargs=waveform_kwargs)

    def select_modes(self, M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, T, **kwargs):
        """
        Selects the harmonic modes that are used for calculating the waveform

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param r0: initial orbital separation of the two objects
        :type r0: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param T: Duration of the waveform in years
        :type T: double, optional

        :rtype: 1d-array[tuples(doubles)]
        """
        return self.waveform_generator.select_modes(M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, T, *args, **kwargs)

    def __call__(self, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, df = None, fmax = None, frequencies = None, **kwargs):
        """
        Calculate the Fourier transform of the plus and cross polarizations of the gravitational wave strain

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param r0: initial orbital separation of the two objects
        :type r0: double
        :param dist: luminosity distance to the source in Gpc
        :type dist: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param dt: Spacing of time samples for corresponding time-domain waveform in seconds. Default is 10. This option is only considered if df, fmax, and frequencies are None.
        :type dt: double, optional
        :param T: Duration of the observed waveform in years. Default is 1. If None, T is set to be 1/df.
        :type T: double, optional
        :param df: Spacing of the frequency samples in Hertz. Default is None. If df is None, then one must specify values for frequencies or T.
        :type df: double, optional
        :param fmax: Maximum sampled frequency in Hertz. Default is None. If fmax is None, then one must specify values for frequencies or dt.
        :type fmax: double, optional
        :param frequencies: Array of frequency samples in Hertz. Default is None. This option overrides df and fmax.
        :type frequencies: list[double] or ndarray[double], optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        :param include_negative_m: True returns the sum of the positive and negative m-modes for each mode in select_modes
        :type include_negative_m: bool, optional
        :param eps: The tolerance to include modes that are subdominant to the power in the (2,2)-mode.
        :type max_samples: double, optional
        :param eps: The tolerance to include modes that are subdominant to the power in the (2,2)-mode.
        :type eps: double, optional

        :rtype: list[two 1d-arrays[double]]
        """
        if "T" in kwargs.keys():
            T = kwargs["T"]
        if "dt" in kwargs.keys():
            dt = kwargs["dt"]
        if "df" in kwargs.keys():
            df = kwargs["df"]
        if "fmax" in kwargs.keys():
            fmax = kwargs["fmax"]
        if "frequencies" in kwargs.keys():
            frequencies = kwargs["frequencies"]

        frequencies, dt, T = sort_frequency_and_time_sampling_arguments(df, fmax, frequencies, dt, T)

        include_negative_m = True
        if "include_negative_m" in kwargs.keys():
            include_negative_m = kwargs["include_negative_m"]

        if "select_modes" in kwargs.keys():
            lmodes = []
            mmodes = []
            if include_negative_m:
                for mode in kwargs["select_modes"]:
                    if mode[1] > 0: # if include_negative_m is True then only keep positive m
                        lmodes.append(mode[0])
                        mmodes.append(mode[1])
                    else:
                        warnings.warn("Warning: Only keeping modes in select_modes with m > 0. Set include_negative_m = False to keep m < 0 modes.")
            else: # if include_negative_m is False then keep all m
                for mode in kwargs["select_modes"]:
                    lmodes.append(mode[0])
                    mmodes.append(mode[1])
            l = np.ascontiguousarray(lmodes, dtype=np.intc)
            m = np.ascontiguousarray(mmodes, dtype=np.intc)
            h = self.waveform_generator.waveform_harmonics(l, m, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, frequencies, dt, T, *args, **kwargs)
        else:
            h = self.waveform_generator.waveform(M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, frequencies, dt, T, *args, **kwargs)
        return h
    
    def harmonics(self, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, df = None, fmax = None, frequencies = None, **kwargs):
        """
        Calculate the spin-weighted spherical harmonic (l, m)-modes of the
        Fourier transform of the plus and cross polarizations of the gravitational wave strain

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param r0: initial orbital separation of the two objects
        :type r0: double
        :param dist: luminosity distance to the source in Gpc
        :type dist: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param dt: Spacing of time samples for corresponding time-domain waveform in seconds. Default is 10. This option is only considered if df, fmax, and frequencies are None.
        :type dt: double, optional
        :param T: Duration of the observed waveform in years. Default is 1. If None, T is set to be 1/df.
        :type T: double, optional
        :param df: Spacing of the frequency samples in Hertz. Default is None. If df is None, then one must specify values for frequencies or T.
        :type df: double, optional
        :param fmax: Maximum sampled frequency in Hertz. Default is None. If fmax is None, then one must specify values for frequencies or dt.
        :type fmax: double, optional
        :param frequencies: Array of frequency samples in Hertz. Default is None. This option overrides df and fmax.
        :type frequencies: list[double] or ndarray[double], optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        :param include_negative_m: True returns the sum of the positive and negative m-modes for each mode in select_modes
        :type include_negative_m: bool, optional
        
        :rtype: list[two 2d-arrays[double]]

        """
        if "T" in kwargs.keys():
            T = kwargs["T"]
        if "dt" in kwargs.keys():
            dt = kwargs["dt"]
        if "df" in kwargs.keys():
            df = kwargs["df"]
        if "fmax" in kwargs.keys():
            fmax = kwargs["fmax"]
        if "frequencies" in kwargs.keys():
            frequencies = kwargs["frequencies"]

        frequencies, dt, T = sort_frequency_and_time_sampling_arguments(df, fmax, frequencies, dt, T)

        include_negative_m = True
        if "include_negative_m" in kwargs.keys():
            include_negative_m = kwargs["include_negative_m"]

        if "select_modes" in kwargs.keys():
            lmodes = []
            mmodes = []
            if include_negative_m:
                for mode in kwargs["select_modes"]:
                    if mode[1] > 0: # if include_negative_m is True then only keep positive m
                        lmodes.append(mode[0])
                        mmodes.append(mode[1])
                    else:
                        warnings.warn("Warning: Only keeping modes in select_modes with m > 0. Set include_negative_m = False to keep m < 0 modes.")
            else: # if include_negative_m is False then keep all m
                for mode in kwargs["select_modes"]:
                    lmodes.append(mode[0])
                    mmodes.append(mode[1])
            l = np.ascontiguousarray(lmodes, dtype=np.intc)
            m = np.ascontiguousarray(mmodes, dtype=np.intc)
            h = self.waveform_generator.waveform_harmonics_grid(l, m, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, frequencies, dt, T, *args, **kwargs)
        else:
            h = self.waveform_generator.waveform_select_harmonics_grid(M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, frequencies, dt, T, *args, **kwargs)
        return h
    
    def harmonics_data(self, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, df = None, fmax = None, frequencies = None, **kwargs):
        """
        Calculate the phases and amplitudes for the spin-weighted spherical harmonic (l, m)-modes of the
        Fourier transform of the plus and cross polarizations of the gravitational wave strain. The
        function returns mode data for both positive and negative m-modes. Given a list of modes [(l_1, m_1), (l_2, m_2), ... , (l_N, m_N)]
        the function returns data in the order [(l_1, m_1), (l_1, -m_1), (l_2, m_2), (l_2, -m_2), ..., (l_N, m_N), (l_N, -m_N)].

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param r0: initial orbital separation of the two objects
        :type r0: double
        :param dist: luminosity distance to the source in Gpc
        :type dist: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param dt: Spacing of time samples for corresponding time-domain waveform in seconds. Default is 10. This option is only considered if df, fmax, and frequencies are None.
        :type dt: double, optional
        :param T: Duration of the observed waveform in years. Default is 1. If None, T is set to be 1/df.
        :type T: double, optional
        :param df: Spacing of the frequency samples in Hertz. Default is None. If df is None, then one must specify values for frequencies or T.
        :type df: double, optional
        :param fmax: Maximum sampled frequency in Hertz. Default is None. If fmax is None, then one must specify values for frequencies or dt.
        :type fmax: double, optional
        :param frequencies: Array of frequency samples in Hertz. Default is None. This option overrides df and fmax.
        :type frequencies: list[double] or ndarray[double], optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        :param include_negative_m: True returns the sum of the positive and negative m-modes for each mode in select_modes
        :type include_negative_m: bool, optional
        
        :rtype: list[two 2d-arrays[double]]

        """
        if "T" in kwargs.keys():
            T = kwargs["T"]
        if "dt" in kwargs.keys():
            dt = kwargs["dt"]
        if "df" in kwargs.keys():
            df = kwargs["df"]
        if "fmax" in kwargs.keys():
            fmax = kwargs["fmax"]
        if "frequencies" in kwargs.keys():
            frequencies = kwargs["frequencies"]

        frequencies, dt, T = sort_frequency_and_time_sampling_arguments(df, fmax, frequencies, dt, T)

        if "select_modes" in kwargs.keys():
            lmodes = []
            mmodes = []
            for mode in kwargs["select_modes"]:
                if mode[1] > 0:
                    lmodes.append(mode[0])
                    mmodes.append(mode[1])
                else:
                    warnings.warn("Warning: Only keeping modes in select_modes with m > 0.")
            l = np.ascontiguousarray(lmodes, dtype=np.intc)
            m = np.ascontiguousarray(mmodes, dtype=np.intc)
            h = self.waveform_generator.waveform_harmonics_phase_amplitude(l, m, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, frequencies, dt, T, *args, **kwargs)
        else:
            h = self.waveform_generator.waveform_select_harmonics_phase_amplitude(M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, frequencies, dt, T, *args, **kwargs)
        return h
    
    def source_frame(self, M, mu, a, r0, theta, phi, Phi_phi0, dt, T, df = None, fmax = None, frequencies = None, **kwargs):
        """
        Calculate the Fourier transform of the plus and cross polarizations of the
        scaled gravitational wave strain :math:`r\\times h/\\mu` in the source frame

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param r0: initial radial separation
        :type r0: double
        :param theta: polar angle of the observor with respect to the Kerr spin vector
        :type theta: double
        :param phi: azimuthal angle of the observor with respect to the Kerr spin vector
        :type phi: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param dt: Spacing of time samples for corresponding time-domain waveform in seconds. Default is 10. If None, dt is set to be 1/(2 fmax).
        :type dt: double, optional
        :param T: Duration of the observed waveform in years. Default is 1. If None, T is set to be 1/df.
        :type T: double, optional
        :param df: Spacing of the frequency samples in Hertz. Default is None. If df is None, then one must specify values for frequencies or T.
        :type df: double, optional
        :param fmax: Maximum sampled frequency in Hertz. Default is None. If fmax is None, then one must specify values for frequencies or dt.
        :type fmax: double, optional
        :param frequencies: Array of frequency samples in Hertz. Default is None. This option overrides df and fmax.
        :type frequencies: list[double] or ndarray[double], optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        
        :rtype: list[two 1d-arrays[double]]

        """
        if "T" in kwargs.keys():
            T = kwargs["T"]
        if "dt" in kwargs.keys():
            dt = kwargs["dt"]
        if "df" in kwargs.keys():
            df = kwargs["df"]
        if "fmax" in kwargs.keys():
            fmax = kwargs["fmax"]
        if "frequencies" in kwargs.keys():
            frequencies = kwargs["frequencies"]

        frequencies, dt, T = sort_frequency_and_time_sampling_arguments(df, fmax, frequencies, dt, T)

        if "select_modes" in kwargs.keys():
            lmodes = []
            mmodes = []
            for mode in kwargs["select_modes"]:
                lmodes.append(mode[0])
                mmodes.append(mode[1])
            l = np.ascontiguousarray(lmodes, dtype=np.intc)
            m = np.ascontiguousarray(mmodes, dtype=np.intc)
            h = self.waveform_generator.waveform_harmonics_source_frame(l, m, M, mu, a, r0, theta, phi, Phi_phi0, frequencies, dt, T, *args, **kwargs)
        else:
            h = self.waveform_generator.waveform_source_frame(M, mu, a, r0, theta, phi, Phi_phi0, frequencies, dt, T, *args, **kwargs)
        return h
    
def sort_frequency_and_time_sampling_arguments(df, fmax, frequencies, dt, T):
    if df is None and frequencies is None and T is None:
            raise ValueError("Values for df, frequencies, and T cannot all be None.")
    
    if fmax is None and frequencies is None and dt is None:
        raise ValueError("Values for fmax, frequencies, and dt cannot all be None.")
    
    # dt is only relevant if it is set to a numerical value and none of the frequency data is specified
    if dt is None or fmax is not None or df is not None or frequencies is not None:
        dt = 0.
    
    if frequencies is not None:
        if not isinstance(frequencies, list) and not isinstance(frequencies, np.ndarray):
            raise ValueError("frequencies must be a list or a ndarray.")
        else:
            frequencies = np.ascontiguousarray(frequencies)
        
        df = np.min(np.diff(frequencies))
        fmax = np.max(frequencies)

        if T is None:
            T = 1/df
    
    if frequencies is None:
        if df is None:
            df = 1/years_to_seconds(T)
        if fmax is None:
            fmax = 1/(2*dt)
        if T is None:
            T = seconds_to_years(1/df)

        frequencies = np.arange(0, fmax, df)

    return frequencies, dt, T

class KerrFrequencyWaveform(KerrCircularFrequencyWaveform):
    """
    Class that generates the gravitational waveform produced by an extreme-mass-ratio inspiral in the
    frequency domain using the adiabatic approximation from black hole perturbation theory and the self-force formalism.
    By default, the waveform is generated in the solar system barycenter frame.
    
    Waveform generation is limited to quasi-circular inspirals in Kerr spacetime, but the generator mirrors the generic
    parametrization used in other EMRI waveform generators (e.g., https://bhptoolkit.org/FastEMRIWaveforms/html/user/main.html)

    :param trajectory_data: a TrajectoryData class which holds interpolants of the relevant trajectory data
    :type trajectory_data: TrajectoryData or None, optional
    :param harmonic_data: a HarmonicAmplitudes class which holds interpolants of the harmonic mode amplitudes
    :type harmonic_data: HarmonicAmplitudes or None, optional
    :param num_threads: the number of threads used to evaluate the waveform
    :type num_threads: int or None, optional  
    """
    def __call__(self, M, mu, a, p0, e0, x0, dist, qS, phiS, qK, phiK, Phi_phi0, Phi_r0, Phi_theta0, dt, T, df = None, fmax = None, frequencies = None, **kwargs):
        """
        Calculate the Fourier transform of the plus and cross polarizations of the gravitational wave strain

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param p0: initial semi-latus rectum
        :type p0: double
        :param e0: initial orbital eccentricity
        :type e0: double
        :param x0: intial cosine of the orbital inclination
        :type x0: double
        :param dist: luminosity distance to the source in Gpc
        :type dist: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param Phi_r0: Phase describing the initial radial position and velocity of the small compact object
        :type Phi_r0: double
        :param Phi_theta0: Phase describing the initial polar position and velocity of the small compact object
        :type Phi_theta0: double
        :param dt: Spacing of time samples for corresponding time-domain waveform in seconds. Default is 10. If None, dt is set to be 1/(2 fmax).
        :type dt: double, optional
        :param T: Duration of the observed waveform in years. Default is 1. If None, T is set to be 1/df.
        :type T: double, optional
        :param df: Spacing of the frequency samples in Hertz. Default is None. If df is None, then one must specify values for frequencies or T.
        :type df: double, optional
        :param fmax: Maximum sampled frequency in Hertz. Default is None. If fmax is None, then one must specify values for frequencies or dt.
        :type fmax: double, optional
        :param frequencies: Array of frequency samples in Hertz. Default is None. This option overrides df and fmax.
        :type frequencies: list[double] or ndarray[double], optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        :param include_negative_m: True returns the sum of the positive and negative m-modes for each mode in select_modes
        :type include_negative_m: bool, optional
        
        :rtype: list[two 1d-arrays[double]]

        """
        a, x0, qK, phiK, Phi_phi0 = check_ssb_frame_parameters(a, x0, qK, phiK, Phi_phi0)
        return super().__call__(M, mu, a, p0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, df = df, fmax = fmax, frequencies = frequencies, **kwargs)
    
    def harmonics(self, M, mu, a, p0, e0, x0, dist, qS, phiS, qK, phiK, Phi_phi0, Phi_r0, Phi_theta0, dt, T, df = None, fmax = None, frequencies = None, **kwargs):
        """
        Calculate the spin-weighted spherical harmonic (l, m)-modes of the
        Fourier transform of the plus and cross polarizations of the gravitational wave strain

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param p0: initial semi-latus rectum
        :type p0: double
        :param e0: initial orbital eccentricity
        :type e0: double
        :param x0: intial cosine of the orbital inclination
        :type x0: double
        :param dist: luminosity distance to the source in Gpc
        :type dist: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param Phi_r0: Phase describing the initial radial position and velocity of the small compact object
        :type Phi_r0: double
        :param Phi_theta0: Phase describing the initial polar position and velocity of the small compact object
        :type Phi_theta0: double
        :param dt: Spacing of time samples for corresponding time-domain waveform in seconds. Default is 10. If None, dt is set to be 1/(2 fmax).
        :type dt: double, optional
        :param T: Duration of the observed waveform in years. Default is 1. If None, T is set to be 1/df.
        :type T: double, optional
        :param df: Spacing of the frequency samples in Hertz. Default is None. If df is None, then one must specify values for frequencies or T.
        :type df: double, optional
        :param fmax: Maximum sampled frequency in Hertz. Default is None. If fmax is None, then one must specify values for frequencies or dt.
        :type fmax: double, optional
        :param frequencies: Array of frequency samples in Hertz. Default is None. This option overrides df and fmax.
        :type frequencies: list[double] or ndarray[double], optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        :param include_negative_m: True returns the sum of the positive and negative m-modes for each mode in select_modes
        :type include_negative_m: bool, optional
        
        :rtype: list[two 2d-arrays[double]]

        """
        a, x0, qK, phiK, Phi_phi0 = check_ssb_frame_parameters(a, x0, qK, phiK, Phi_phi0)
        return super().harmonics(M, mu, a, p0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, df = df, fmax = fmax, frequencies = frequencies, **kwargs)
    
    def harmonics_data(self, M, mu, a, p0, e0, x0, dist, qS, phiS, qK, phiK, Phi_phi0, Phi_r0, Phi_theta0, dt, T, df = None, fmax = None, frequencies = None, **kwargs):
        """
        Calculate the phases and amplitudes for spin-weighted spherical harmonic (l, m)-modes of the
        Fourier transform of the plus and cross polarizations of the gravitational wave strain

        :param M: mass (in solar masses) of the massive black hole
        :type M: double
        :param mu: mass (in solar masses) of the (smaller) stellar-mass compact object
        :type mu: double
        :param a: dimensionless black hole spin
        :type a: double
        :param p0: initial semi-latus rectum
        :type p0: double
        :param e0: initial orbital eccentricity
        :type e0: double
        :param x0: intial cosine of the orbital inclination
        :type x0: double
        :param dist: luminosity distance to the source in Gpc
        :type dist: double
        :param qS: polar angle of the source's sky location
        :type qS: double
        :param phiS: azimuthal angle of the source's sky location
        :type phiS: double
        :param qK: polar angle of the Kerr spin vector
        :type qK: double
        :param phiK: azimuthal angle of the Kerr spin vector
        :type phiK: double
        :param Phi_phi0: Initial azimuthal position of the small compact object
        :type Phi_phi0: double
        :param Phi_r0: Phase describing the initial radial position and velocity of the small compact object
        :type Phi_r0: double
        :param Phi_theta0: Phase describing the initial polar position and velocity of the small compact object
        :type Phi_theta0: double
        :param dt: Spacing of time samples for corresponding time-domain waveform in seconds. Default is 10. If None, dt is set to be 1/(2 fmax).
        :type dt: double, optional
        :param T: Duration of the observed waveform in years. Default is 1. If None, T is set to be 1/df.
        :type T: double, optional
        :param df: Spacing of the frequency samples in Hertz. Default is None. If df is None, then one must specify values for frequencies or T.
        :type df: double, optional
        :param fmax: Maximum sampled frequency in Hertz. Default is None. If fmax is None, then one must specify values for frequencies or dt.
        :type fmax: double, optional
        :param frequencies: Array of frequency samples in Hertz. Default is None. This option overrides df and fmax.
        :type frequencies: list[double] or ndarray[double], optional

        :param pad_output: True returns the waveform for the full duration T years even if the system merges before T years has elasped
        :type pad_output: bool, optional
        :param select_modes: A list of tuples :math:`(l, m)` that select which modes to include in the waveform calculation
        :type select_modes: list[tuple(double)] or ndarray[tuple(double)], optional
        :param return_list: True returns the plus and cross polarizations of the waveform as separate ndarrays
        :type return_list: bool, optional
        :param include_negative_m: True returns the sum of the positive and negative m-modes for each mode in select_modes
        :type include_negative_m: bool, optional
        
        :rtype: list[two 2d-arrays[double]]

        """
        a, x0, qK, phiK, Phi_phi0 = check_ssb_frame_parameters(a, x0, qK, phiK, Phi_phi0)
        return super().harmonics_data(M, mu, a, p0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, df = df, fmax = fmax, frequencies = frequencies, **kwargs)

def source_angles(qS, phiS, qK, phiK):
    """
    Calculate the sky location :math:`(\\theta, \\phi)` of the observor in the
    source frame using the sky location and orientation of the source in
    the SSB frame

    :param qS: polar angle of the source's sky location
    :type qS: double
    :param phiS: azimuthal angle of the source's sky location
    :type phiS: double
    :param qK: polar angle of the Kerr spin vector
    :type qK: double
    :param phiK: azimuthal angle of the Kerr spin vector
    :type phiK: double

    :rtype: tuple(double, double)
    """
    phi = -0.5*np.pi
    theta = np.arccos(-(np.sin(qS)*np.sin(qK)*np.cos(phiS - phiK) + np.cos(qS)*np.cos(qK)))

    return (theta, phi)

def polarization(qS, phiS, qK, phiK):
    """
    Calculate the rotation of polarization angle :math:`e^{i\\psi}` due to transforming from
    the plus and cross polarizations in the source frame to the plus and
    cross polarization in the SSB frame.

    :param qS: polar angle of the source's sky location
    :type qS: double
    :param phiS: azimuthal angle of the source's sky location
    :type phiS: double
    :param qK: polar angle of the Kerr spin vector
    :type qK: double
    :param phiK: azimuthal angle of the Kerr spin vector
    :type phiK: double

    :rtype: complex

    """
    real_part = np.cos(qS)*np.sin(qK)*np.cos(phiS - phiK) - np.cos(qK)*np.sin(qS)
    imag_part = -np.sin(qK)*np.sin(phiS - phiK)
    if abs(real_part) + abs(imag_part) == 0.:
        return 1.

    return (real_part + 1.j*imag_part)/(real_part - 1.j*imag_part)

def solar_mass_to_seconds(mass):
    """
    Converts units of solar mass in :math:`G=c=1` convention to units of seconds in the MKS convention
    
    :param mass: mass in units of solar masses
    :type mass: double
    
    :rtype: double
    """
    return mass*Modot_GC1_to_S

def seconds_to_solar_mass(time):
    """
    Converts units of seconds in the MKS convention to units of solar mass in :math:`G=c=1` convention
    
    :param time: time in units of seconds
    :type time: double
    
    :rtype: double
    """
    return time/Modot_GC1_to_S

def solar_mass_to_meters(mass):
    """
    Converts units of solar mass in :math:`G=c=1` convention to units of meters in the MKS convention
    
    :param mass: mass in units of solar masses
    :type mass: double
    
    :rtype: double
    """
    return mass*Modot_GC1_to_M

def solar_mass_to_parsecs(mass):
    """
    Converts units of solar mass in :math:`G=c=1` convention to units of parsecs in the MKS convention
    
    :param mass: mass in units of solar masses
    :type mass: double
    
    :rtype: double
    """
    return mass*Modot_GC1_to_PC

def parsecs_to_solar_mass(length):
    """
    Converts units of parsecs in the MKS convention to units of solar mass in :math:`G=c=1` convention
    
    :param length: length in units of solar masses
    :type length: double
    
    :rtype: double
    """
    return length/Modot_GC1_to_PC

def seconds_to_years(time):
    """
    Converts units of seconds to units of sidereal years
    
    :param time: time in seconds
    :type time: double
    
    :rtype: double
    """
    return time/yr_MKS

def years_to_seconds(time):
    """
    Converts units of sidereal years to units of seconds
    
    :param time: time in sidereal years
    :type time: double
    
    :rtype: double
    """
    return time*yr_MKS

def scaled_amplitude(mu, dist):
    """
    Gives the overall scaling of the waveform amplitude :math:`\\mu/D_L` for a binary with 
    small mass :math:`\\mu` and an observer at distance :math:`D_L`
    
    :param mu: mass of the small body in solar masses
    :type mu: double
    :param dist: distance between the binary and observor in Gpc
    :type dist: double
    
    :rtype: double
    """
    return Modot_GC1_to_PC*mu/(dist*1.e9)

def frequencies(dt, T):
    """
    Frequency sampling of the frequency-domain gravitational wave signal for a given
    time step and signal duration

    :param dt: time step in seconds
    :type dt: double
    :param T: signal duration in years
    :type T: double

    :rtype: array[double]
    """
    samples = int(T*yr_MKS/dt + 1)
    return np.fft.rfftfreq(samples, d=dt)


## Custom Classes ##