#cython: profile=True

from libcpp.vector cimport vector
# from libcpp.complex cimport complex
import numpy as np
cimport numpy as np
from libc.string cimport memcpy
from libc.math cimport sin, cos, acos, abs, M_PI
from cython.operator import dereference
cimport openmp
import os

# G_const = 6.67430e-11 # CODATA value from NIST
# # this is the gravitational parameter which is measured more accurately than Modot or G 
# GM_const = 1.32712440041279419e+20 # new value from https://iopscience.iop.org/article/10.3847/1538-3881/abd414
# Modot_const = 1.98841e+30 # GM_const/G_const and keeping 6 sigfigs
# c_const = 299792458.
# pc_const = 3.0856775814913674e+16
# yr_const = 31558149.763545603 # this is the sidereal year (365.256363004 days)

cdef extern from "waveform.hpp":
    double solar_mass_to_seconds(double)
    double scale_strain_amplitude(double, double)
    int inspiral_time_steps(double, double, double, double, double, double)
    void output_waveform_td(float*, float*, double, double, double, double, double, double, double, double, double, int)
    void output_waveform_harmonic_td(float*, float*, int, int, double, double, double, double, double, double, double, double, double, int)
    void output_trajectory_td(double *, double *, double, double, double, double, double, double, double)
    void output_downsampled_trajectory(double *, double *, double *, int, double, double, double, double, double, double, double, double)
    void output_flux(double *, double *, double *, int)
    double G_const
    double GM_const
    double c_const
    double Modot_const
    double pc_const
    double yr_const

cdef extern from "trajectory.hpp":
    double kerr_isco_radius(double, int)
    double kerr_isco_frequency(double a)
    double kerr_geo_radius_circ(double a, double Omega)
    double kerr_geo_azimuthal_frequency_circ_time(double a, double r, int sgnX)
    double kerr_geo_azimuthal_frequency_circ_time(double a, double r)
    double alpha_of_a_omega(const double &a, const double &omega, const double &oISCO)
    double alpha_of_a_omega(const double &a, const double &omega)
    double dalpha_domega_of_a_omega(const double &a, const double &omega, const double &oISCO)
    double omega_of_a_alpha(const double &a, const double &alpha, const double &oISCO)
    double omega_of_a_alpha(const double &a, const double &alpha)
    double spin_of_chi(const double &chi)
    double chi_of_spin(const double &a)

Modot_SI = Modot_const
GM_SI = GM_const
G_SI = G_const
c_SI = c_const
pc_SI = pc_const
yr_SI = yr_const
OMEGA_MIN = 2.e-3

def solar_mass_to_seconds(double mass):
    return mass*GM_const/c_const**3

def seconds_to_solar_mass(double seconds):
    return seconds/solar_mass_to_seconds(1.)

def solar_mass_to_meters(double mass):
    return mass*GM_const/c_const**2

def solar_mass_to_parsecs(double mass):
    return solar_mass_to_meters(mass)/pc_const

def parsecs_to_solar_mass(double pc):
    return pc/solar_mass_to_parsecs(1.)

def seconds_to_years(double seconds):
    return seconds/yr_const

def years_to_seconds(double years):
    return years*yr_const

def scaled_amplitude(double mu, double dist):
    return scale_strain_amplitude(mu, dist)

def omega_of_a_alpha(a, alpha):
    oISCO13 = kerr_isco_frequency(a)**(1./3.)
    omin13 = OMEGA_MIN**(1./3.)
    return (oISCO13 - alpha**2*(oISCO13 - omin13))**3.

def min_orbital_frequency(a):
    return omega_of_a_alpha(a, 1.)

def max_orbital_frequency(a):
    return omega_of_a_alpha(a, 0.)

def circular_energy(a, r):
    v = 1./np.sqrt(r)
    return (1. - 2.*v**2 + a*v**3)/np.sqrt(1. - 3.*v**2 + 2.*a*v**3)

def circular_frequency(a, r):
    return r**(-1.5)/(1 + a*r**(-1.5))

def kerr_radius(a, omega):
    return ((1. - a*omega)**2/omega**2)**(1./3.)

def kerr_isco(a):
    z1 = 1 + (1 - a**2)**(1./3.)*((1 + a)**(1./3.) + (1 - a)**(1./3.))
    z2 = np.sqrt(3*a**2 + z1**2)

    return 3. + z2 - np.sign(a)*np.sqrt((3. - z1)*(3. + z1 + 2.*z2))

def kerr_isco_frequency(a):
    return circular_frequency(a, kerr_isco(a))

def energy_flux(a, r):
    return energy_flux_of_a_omega(a, circular_frequency(a, r))

def energy_flux_of_a_omega(a, omega):
    cdef int ashape, oshape, nsize
    if isinstance(a, np.ndarray):
        ashape = a.shape[0]
    elif np.asarray(a).shape is ():
        a = np.array([a])
        ashape = 1
    elif isinstance(a, list):
        return energy_flux(np.array(a), omega)

    if isinstance(omega, np.ndarray):
        oshape = omega.shape[0]
    elif np.asarray(omega).shape is ():
        omega = np.array([omega])
        oshape = 1
    elif isinstance(omega, list):
        return energy_flux(a, np.array(omega))

    if ashape > 1 and oshape > 1:
        assert oshape == ashape, "Arguments have unequal shapes of ({},) and ({},)".format(ashape, oshape)

    nsize = max(oshape, ashape)
    cdef np.ndarray[ndim=1, dtype=np.float64_t] dedt = np.empty((nsize), dtype=np.float64)
    
    if ashape > 1 and oshape == 1:
        omega = np.full((nsize), omega[0])
    elif ashape == 1 and oshape > 1:
        a = np.full((nsize), a[0])

    omin = min_orbital_frequency(a)
    omax = kerr_isco_frequency(a)
    if not np.all(omega >= omin):
        raise ValueError("Radii greater than {} are not supported".format(np.max(kerr_radius(a, omin))))
    if not np.all(omega <= omax):
        raise ValueError("Radii less than {} are not supported".format(np.min(kerr_radius(a, omax))))

    return energy_flux_numpy(a, omega)

cdef energy_flux_numpy(np.ndarray[ndim=1, dtype=np.float64_t] a, np.ndarray[ndim=1, dtype=np.float64_t] omega):
    cdef int nsize = a.shape[0]
    cdef np.ndarray[ndim=1, dtype=np.float64_t] dedt = np.empty((nsize), dtype=np.float64)
    cdef np.ndarray[ndim=1, dtype=np.float64_t] ac = a
    cdef np.ndarray[ndim=1, dtype=np.float64_t] omegac = omega

    output_flux(&dedt[0], &ac[0], &omegac[0], nsize)

    return dedt

def downsampled_trajectory(double M, double mu, double a, double r0, double dist, double theta, double phi, double T, double dt, int nsamples=100):
    cdef np.ndarray[ndim=1, dtype=np.float64_t] t = np.zeros(nsamples, dtype=np.float64)
    cdef np.ndarray[ndim=1, dtype=np.float64_t] alpha = np.zeros(nsamples, dtype=np.float64)
    cdef np.ndarray[ndim=1, dtype=np.float64_t] phase = np.zeros(nsamples, dtype=np.float64)

    output_downsampled_trajectory(&t[0], &alpha[0], &phase[0], nsamples, theta, phi, mu, M, a, r0, T, dt)
    return (t, alpha, phase)

def trajectory(double M, double mu, double a, double r0, double Phi_phi0, double T, double dt):
    cdef int time_steps 
    time_steps = inspiral_time_steps(mu, M, a, r0, T, dt)
    cdef np.ndarray[ndim=1, dtype=np.float64_t] rp = np.zeros(time_steps, dtype=np.float64)
    cdef np.ndarray[ndim=1, dtype=np.float64_t] phase = np.zeros(time_steps, dtype=np.float64)

    output_trajectory_td(&rp[0], &phase[0], mu, M, a, r0, Phi_phi0, T, dt)

    return (rp, phase)

def waveform(double M, double mu, double a, double r0, double dist, double theta, double phi, double T, double dt, int num_threads=0, int pad_output=0):
    cdef int time_steps 
    if pad_output:
        time_steps = int(T*yr_const/dt) + 1
    else:
        time_steps = inspiral_time_steps(mu, M, a, r0, T, dt)
    cdef np.ndarray[ndim=1, dtype=np.float32_t] hp = np.zeros(time_steps, dtype=np.float32)
    cdef np.ndarray[ndim=1, dtype=np.float32_t] hc = np.zeros(time_steps, dtype=np.float32)
    
    output_waveform_td(&hp[0], &hc[0], mu, M, a, r0, T, dist, theta, phi, dt, num_threads)

    return hp - 1j*hc

def waveform_harmonic(int l, int m, double M, double mu, double a, double r0, double dist, double theta, double phi, double T, double dt, int num_threads=0):
    cdef int time_steps 
    time_steps = inspiral_time_steps(mu, M, a, r0, T, dt)
    cdef np.ndarray[ndim=1, dtype=np.float32_t] hp = np.zeros(time_steps, dtype=np.float32)
    cdef np.ndarray[ndim=1, dtype=np.float32_t] hc = np.zeros(time_steps, dtype=np.float32)
    
    output_waveform_harmonic_td(&hp[0], &hc[0], l, m, mu, M, a, r0, T, dist, theta, phi, dt, num_threads)

    return hp - 1j*hc

def source_angles(double qk, double phik, double qs, double phis):
    cdef double phi = -0.5*M_PI
    cdef double theta = acos(-(sin(qs)*sin(qk)*cos(phis - phik) + cos(qs)*cos(qk)))

    return (theta, phi)

def polarization(double qk, double phik, double qs, double phis):
    cdef double real_part = cos(qs)*sin(qk)*cos(phis - phik) - cos(qk)*sin(qs)
    cdef double imag_part = -sin(qk)*sin(phis - phik)
    if abs(real_part) + abs(imag_part) == 0.:
        return 0.j

    return (real_part + 1.j*imag_part)**2/(real_part**2 + imag_part**2)

def get_max_threads():
    return openmp.omp_get_max_threads()

def get_num_threads():
    return openmp.omp_get_max_threads()

def alphaOfSpinOmega(double a, double omega):
    return alpha_of_a_omega(a, omega)

def omegaOfSpinAlpha(double a, double alpha):
    return omega_of_a_alpha(a, alpha)

class KerrCircularWaveform:
    def __init__(self, nthreads=None, pad_output=False, convert_float=False):
        if nthreads is None:
            self.nthreads=0
        else:
            self.nthreads=nthreads

        if pad_output:
            self.pad_output = 1
        else:
            self.pad_output = 0
        
        self.convert_float = convert_float

    def __call__(self,
        double M,
        double mu,
        double a,
        double p0,
        double e0,
        double x0,
        double dist,
        double qS,
        double phiS,
        double qK,
        double phiK,
        double Phi_phi0,
        double Phi_theta0,
        double Phi_r0,
        T=1.,
        dt=10.
    ): 
        theta, phi = source_angles(qK, phiK, qS, phiS)
        rot = polarization(qK, phiK, qS, phiS)

        if self.convert_float:
            return rot*(waveform(M, mu, a, p0, dist, theta, phi, T, dt, num_threads=self.nthreads, pad_output=self.pad_output).astype(np.complex128, copy=False))
        else:
            return rot*waveform(M, mu, a, p0, dist, theta, phi, T, dt, num_threads=self.nthreads, pad_output=self.pad_output)



