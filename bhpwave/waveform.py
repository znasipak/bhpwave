import numpy as np
from bhpwaveformcy import (WaveformHarmonicGeneratorPyWrapper,
                           WaveformFourierGeneratorPy,
                           WaveformGeneratorPy,
                           TrajectoryDataPy,
                           InspiralGeneratorPy,
                           HarmonicAmplitudesPy)

import os
from bhpwave.constants import *

path_to_file = os.path.dirname(os.path.abspath(__file__))
traj_path = path_to_file + "/data/trajectory.txt"
amplitude_path = path_to_file + "/data/circ_data"

try:
    CPU_MAX = len(os.sched_getaffinity(0))
except:
    CPU_MAX = os.cpu_count()

class KerrCircularWaveformBase:
    """
    Base class that generates a gravitational waveform produced by an extreme-mass-ratio inspiral
    using the adiabatic approximation from black hole perturbation theory and the self-force formalism.
    The waveform is generated in units :math:`G = c = 1`, with time measured in units of :math:`M`,
    where :math:`M` is the mass of the more massive body. 

    Waveform generation is limited to quasi-circular inspirals in Kerr spacetime.
    
    :param trajectory_data: a TrajectoryDataPy class which holds interpolants of the relevant trajectory data
    :type trajectory_data: TrajectoryDataPy or None, optional
    :param harmonic_data: a HarmonicAmplitudesPy class which holds interpolants of the harmonic mode amplitudes
    :type harmonic_data: TrajectoryDataPy or None, optional
    :param num_threads: the number of threads used to evaluate the waveform
    :type num_threads: int or None, optional
    """
    def __init__(self, trajectory_data=None, harmonic_data=None, num_threads=None):
        if num_threads is None:
            num_threads = CPU_MAX
        if trajectory_data is None:
            self.trajectory_data = TrajectoryDataPy(filename=traj_path,dealloc_flag=False)
        else:
            self.trajectory_data = trajectory_data
        if harmonic_data is None:
            self.harmonic_data = HarmonicAmplitudesPy(filebase=amplitude_path,dealloc_flag=False)
        else:
            self.harmonic_data = harmonic_data

        self.inspiral_generator = InspiralGeneratorPy(self.trajectory_data, num_threads=num_threads)
        waveform_kwargs = {
            "num_threads": num_threads
        }
        self.waveform_generator = WaveformHarmonicGeneratorPyWrapper(self.harmonic_data, waveform_kwargs=waveform_kwargs)

    def generate_base_waveform(self, massratio, a, r0, dt, T, theta, phi, **kwargs):
        """
        Calculate the complex gravitational wave strain :math:`h = h_+ - i h_\\times` measured in the
        solar system barycenter (SSB) frame

        :param massratio: dimensionless ratio between the smaller and larger masses of the binary
        :type massratio: double
        :param a: dimensionless spin of the massive black hole (MBH)
        :type a: double
        :param r0: initial orbital separation of the two objects
        :type r0: double
        :param dt: Spacing of time samples in units of mass of the MBH
        :type dt: double
        :param T: Duration of the waveform in units of mass of the MBH
        :type T: double
        :param theta: polar angle of the observor's sky location with respect to the MBH's spin vector
        :type theta: double
        :param phi: azimuthal angle of the observor's sky location
        :type phi: double

        :rtype: 1d-array[complex]
        """
        if "num_threads" in kwargs.keys():
            inspiral = self.inspiral_generator(massratio, a, r0, dt, T, num_threads=kwargs["num_threads"])
        else:
            inspiral = self.inspiral_generator(massratio, a, r0, dt, T)

        if "modes" in kwargs.keys():
            modes = np.array(kwargs["modes"])
            if np.array(kwargs["modes"]).shape[0] == 2:
                h = self.waveform_generator(modes[0], modes[1], inspiral, theta, phi)
            else:
                raise KeyError("modes key must be an array or list with dimension (2 x # of modes)")
        else:
            h = self.waveform_generator(inspiral, theta, phi)
    
        return h.plus - 1.j*h.cross
    
    def __call__(self, massratio, a, r0, dt, T, theta, phi, **kwargs):
        """
        Calculate the complex gravitational wave strain :math:`h = h_+ - i h_\\times` measured in the
        solar system barycenter (SSB) frame.

        :param massratio: dimensionless ratio between the smaller and larger masses of the binary
        :type massratio: double
        :param a: dimensionless spin of the massive black hole (MBH)
        :type a: double
        :param r0: initial orbital separation of the two objects
        :type r0: double
        :param dt: Spacing of time samples in units of mass of the MBH
        :type dt: double
        :param T: Duration of the waveform in units of mass of the MBH
        :type T: double
        :param theta: polar angle of the observor's sky location with respect to the MBH's spin vector
        :type theta: double
        :param phi: azimuthal angle of the observor's sky location
        :type phi: double

        :rtype: 1d-array[complex]
        """
        return self.generate_base_waveform(massratio, a, r0, dt, T, theta, phi, **kwargs)

class KerrCircularWaveform:
    """
    Class that generates the gravitational waveform produced by an extreme-mass-ratio inspiral
    using the adiabatic approximation from black hole perturbation theory and the self-force formalism.
    By default, the waveform is generated in the solar system barycenter frame.

    Waveform generation is limited to quasi-circular inspirals in Kerr spacetime, but the generator mirrors the generic
    parametrization used in other EMRI waveform generators (e.g., https://bhptoolkit.org/FastEMRIWaveforms/html/user/main.html)

    :param trajectory_data: a TrajectoryDataPy class which holds interpolants of the relevant trajectory data
    :type trajectory_data: TrajectoryDataPy or None, optional
    :param harmonic_data: a HarmonicAmplitudesPy class which holds interpolants of the harmonic mode amplitudes
    :type harmonic_data: TrajectoryDataPy or None, optional
    :param num_threads: the number of threads used to evaluate the waveform
    :type num_threads: int or None, optional
    """
    def __init__(self, trajectory_data=None, harmonic_data=None, num_threads=None, frequency_domain=False):
        if num_threads is None:
            num_threads = CPU_MAX
        if trajectory_data is None:
            self.trajectory_data = TrajectoryDataPy(filename=traj_path, dealloc_flag=False)
        else:
            self.trajectory_data = trajectory_data
        if harmonic_data is None:
            self.harmonic_data = HarmonicAmplitudesPy(filebase=amplitude_path, dealloc_flag=False)
        else:
            self.harmonic_data = harmonic_data

        waveform_kwargs = {
            "num_threads": num_threads
        }
        if frequency_domain:
            self.waveform_generator = WaveformFourierGeneratorPy(self.trajectory_data, self.harmonic_data, waveform_kwargs=waveform_kwargs)
        else:
            self.waveform_generator = WaveformGeneratorPy(self.trajectory_data, self.harmonic_data, waveform_kwargs=waveform_kwargs)

    def select_modes(self, M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, dt=10., T=1., **kwargs):
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
        return self.waveform_generator.select_modes(M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, dt, T, **kwargs)

    def __call__(self, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt=10., T=1., **kwargs):
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

        :rtype: 1d-array[complex]
        """
        h = self.waveform_generator.waveform(M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, **kwargs)
        return h.plus - 1.j*h.cross

class KerrWaveform(KerrCircularWaveform):
    """
    Class that generates the gravitational waveform produced by an extreme-mass-ratio inspiral
    using the adiabatic approximation from black hole perturbation theory and the self-force formalism.
    By default, the waveform is generated in the solar system barycenter frame.
    
    Waveform generation is limited to quasi-circular inspirals in Kerr spacetime, but the generator mirrors the generic
    parametrization used in other EMRI waveform generators (e.g., https://bhptoolkit.org/FastEMRIWaveforms/html/user/main.html)

    :param trajectory_data: a TrajectoryDataPy class which holds interpolants of the relevant trajectory data
    :type trajectory_data: TrajectoryDataPy or None, optional
    :param harmonic_data: a HarmonicAmplitudesPy class which holds interpolants of the harmonic mode amplitudes
    :type harmonic_data: TrajectoryDataPy or None, optional
    :param num_threads: the number of threads used to evaluate the waveform
    :type num_threads: int or None, optional  

    """
    def __call__(self, M, mu, a, p0, e0, x0, dist, qS, phiS, qK, phiK, Phi_phi0, Phi_r0, Phi_theta0, dt=10., T=1., **kwargs):
        """
        Calculate the complex gravitational wave strain

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
        
        :rtype: 1d-array[complex] or list[two 1d-arrays[double]]

        """
        if "T" in kwargs.keys():
            T = kwargs["T"]
        if "dt" in kwargs.keys():
            dt = kwargs["dt"]

        if "select_modes" in kwargs.keys():
            lmodes = []
            mmodes = []
            for mode in kwargs["select_modes"]:
                lmodes.append(mode[0])
                mmodes.append(mode[1])
            l = np.ascontiguousarray(lmodes)
            m = np.ascontiguousarray(mmodes)
            h = self.waveform_generator.waveform_harmonics(l, m, M, mu, a, p0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, **kwargs)
        else:
            h = self.waveform_generator.waveform(M, mu, a, p0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, **kwargs)
        return h
    
    def source_frame(self, M, mu, a, p0, theta, phi, Phi_phi0, dt=10., T=1., **kwargs):
        if "T" in kwargs.keys():
            T = kwargs["T"]
        if "dt" in kwargs.keys():
            dt = kwargs["dt"]

        if "select_modes" in kwargs.keys():
            lmodes = []
            mmodes = []
            for mode in kwargs["select_modes"]:
                lmodes.append(mode[0])
                mmodes.append(mode[1])
            l = np.ascontiguousarray(lmodes)
            m = np.ascontiguousarray(mmodes)
            h = self.waveform_generator.waveform_harmonics_source_frame(l, m, M, mu, a, p0, theta, phi, Phi_phi0, dt, T, **kwargs)
        else:
            h = self.waveform_generator.waveform_source_frame(M, mu, a, p0, theta, phi, Phi_phi0, dt, T, **kwargs)
        return h

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
        return 0.j

    return (real_part + 1.j*imag_part)**2/(real_part**2 + imag_part**2)

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

def frequencies(T, dt):
    samples = int(T*yr_MKS/dt + 1)
    return np.fft.rfftfreq(samples, d=dt)