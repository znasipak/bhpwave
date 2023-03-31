import numpy as np
cimport numpy as np
from libcpp.string cimport string
from libcpp.pair cimport pair
from cython.operator import dereference
import os
import warnings

cdef unicode default_harmonic_filebase = '../bhpwave/data/circ_data'

include "trajectory_wrap.pyx"

cdef extern from "spline.hpp":
    cdef cppclass EigenCubicInterpolator:
        pass

    cdef cppclass EigenBicubicInterpolator:
        EigenCubicInterpolator reduce_x(const double x)
        EigenCubicInterpolator reduce_y(const double y)
        pass

cdef extern from "harmonics.hpp":
    cdef cppclass HarmonicSpline:
        HarmonicSpline(double spin, EigenCubicInterpolator amplitude_spline, EigenCubicInterpolator phase_spline)

        double amplitude(double alpha)
        double phase(double alpha)

        double amplitude_of_omega(double omega)
        double phase_of_omega(double omega)
        double phase_of_omega_derivative(double omega)

    cdef cppclass HarmonicSpline2D:
        HarmonicSpline2D(int j, int m, string filebase)

        double amplitude(double chi, double alpha)
        double phase(double chi, double alpha)

        double amplitude_of_a_omega(double a, double omega)
        double phase_of_a_omega(double a, double omega)
        double phase_of_a_omega_derivative(double a, double omega)

    cdef cppclass HarmonicAmplitudes:
        HarmonicAmplitudes(int lmodes[], int mmodes[], int modeNum, string filepath_base)

        double amplitude(int l, int m, double chi, double alpha)
        double phase(int l, int m, double chi, double alpha)
        int key_check(pair[int, int] key)

        double amplitude_of_a_omega(int l, int m, double a, double omega)
        double phase_of_a_omega(int l, int m, double a, double omega)
        double phase_of_a_omega_derivative(int l, int m, double a, double omega)

    cdef cppclass HarmonicOptions:
        HarmonicOptions()
        HarmonicOptions(double eps, int max)
        double epsilon
        int max_samples

    cdef cppclass HarmonicModeContainer:
        HarmonicModeContainer()
        vector[int] lmodes
        vector[int] mmodes
        vector[double] plusY
        vector[double] crossY

    cdef cppclass HarmonicSelector:
        HarmonicSelector(HarmonicAmplitudes &harm, HarmonicOptions opts)

        double modePower(int l, int m, InspiralContainer &inspiral)
        int gradeMode(int l, int m, InspiralContainer &inspiral, double power22)
        int gradeMode(int l, int m, InspiralContainer &inspiral, double power22, double &plusYlm, double &crossYlm, double theta)
        HarmonicModeContainer selectModes(InspiralContainer &inspiral, double theta)

        double modePower(int l, int m, InspiralContainer &inspiral, HarmonicOptions opts)
        int gradeMode(int l, int m, InspiralContainer &inspiral, double power22, HarmonicOptions opts)
        int gradeMode(int l, int m, InspiralContainer &inspiral, double power22, double &plusYlm, double &crossYlm, double theta, HarmonicOptions opts)
        HarmonicModeContainer selectModes(InspiralContainer &inspiral, double theta, HarmonicOptions opts)

cdef class HarmonicAmplitude2D:
    cdef HarmonicSpline2D *harmcpp

    def __init__(self, int j, int m, unicode filebase=default_harmonic_filebase):
        self.harmcpp = new HarmonicSpline2D(j, m, filebase.encode())

    def __dealloc__(self):
        del self.harmcpp

    def amplitude(self, double a, double omega):
        return self.harmcpp.amplitude_of_a_omega(a, omega)

    def phase(self, double a, double omega):
        return self.harmcpp.phase_of_a_omega(a, omega)

    def phase_omega_derivative(self, double a, double omega):
        return self.harmcpp.phase_of_a_omega_derivative(a, omega)

# cdef class HarmonicAmplitude1D:
#     cdef HarmonicSpline *harmcpp

#     def __cinit__(self, double a, HarmonicAmplitude2D harm2D):
#         self.harmcpp = new HarmonicSpline(a, harm2D.harmcpp.getReducedAmplitudeSpline(a), harm2D.harmcpp.getReducedPhaseSpline(a))

#     def __dealloc__(self):
#         del self.harmcpp

#     def amplitude(self, double omega):
#         return self.harmcpp.amplitude_of_omega(omega)

#     def phase(self, double omega):
#         return self.harmcpp.phase_of_omega(omega)

#     def phase_omega_derivative(self, double omega):
#         return self.harmcpp.phase_of_omega_derivative(omega)

cdef int[::1] default_lmodes():
    cdef int lmax = 15
    cdef int modeNum = 0
    for l in range(2, lmax+1):
        modeNum += l
    cdef int[::1] lmodes = np.empty(modeNum, dtype=np.int32)
    cdef int i = 0
    for l in range(2, lmax+1):
        for m in range(1, l+1):
            lmodes[i] = l
            i += 1
    return lmodes

cdef int[::1] default_mmodes():
    cdef int lmax = 15
    cdef int modeNum = 0
    for l in range(2, lmax+1):
        modeNum += l
    cdef int[::1] mmodes = np.empty(modeNum, dtype=np.int32)
    cdef int i = 0
    for l in range(2, lmax+1):
        for m in range(1, l+1):
            mmodes[i] = m
            i += 1
    return mmodes

cdef int[::1] DEFAULT_LMODES = default_lmodes()
cdef int[::1] DEFAULT_MMODES = default_mmodes()

cdef class HarmonicModeContainerWrapper:
    cdef HarmonicModeContainer modecpp

    cdef wrap(self, HarmonicModeContainer mode):
        self.modecpp = mode

    @property
    def lmodes(self):
        cdef int[::1] arr = <int [:self.size]>self.modecpp.lmodes.data()
        return np.asarray(arr)

    @property
    def mmodes(self):
        cdef int[::1] arr = <int [:self.size]>self.modecpp.mmodes.data()
        return np.asarray(arr)
    
    @property
    def modes(self):
        return np.array([self.lmodes, self.mmodes])

    @property
    def size(self):
        return self.modecpp.lmodes.size()

    @property
    def modecount(self):
        return self.size


cdef class HarmonicAmplitudesPy:
    cdef HarmonicAmplitudes *harmonicscpp
    cdef bint dealloc_flag

    def __cinit__(self, int[::1] lmodes = DEFAULT_LMODES, int[::1] mmodes = DEFAULT_MMODES, unicode filebase = default_harmonic_filebase, bint dealloc_flag = True):
        self.harmonicscpp = new HarmonicAmplitudes(&lmodes[0], &mmodes[0], lmodes.shape[0], filebase.encode())
        self.dealloc_flag = dealloc_flag

    def __dealloc__(self):
        if self.dealloc_flag:
            warnings.warn("Deallocating HarmonicAmplitudesPy object", UserWarning)
        del self.harmonicscpp
    
    def amplitude(self, int l, int m, double a, double r):
        return self.harmonicscpp.amplitude_of_a_omega(l, m, a, abs(kerr_geo_orbital_frequency_circ(a, r)))

    def phase(self, int l, int m, double a, double r):
        return self.harmonicscpp.phase_of_a_omega(l, m, a, abs(kerr_geo_orbital_frequency_circ(a, r)))

    def amplitude_array(self, int l, int m, a, r):
        a_array = np.ascontiguousarray(a)
        r_array = np.ascontiguousarray(r)
        if a_array.shape == r_array.shape:
            amp_list = [self.harmonicscpp.amplitude_of_a_omega(l, m, a_array[i], abs(kerr_geo_orbital_frequency_circ(a_array[i], r_array[i]))) for i in range(a_array.shape[0])]
        else:
            amp_list = [[self.harmonicscpp.amplitude_of_a_omega(l, m, aa, abs(kerr_geo_orbital_frequency_circ(aa, rr))) for rr in r_array] for aa in a_array]
        return np.array(amp_list)

    def phase_array(self, int l, int m, a, r):
        a_array = np.ascontiguousarray(a)
        r_array = np.ascontiguousarray(r)
        if a_array.shape == r_array.shape:
            amp_list = [self.harmonicscpp.phase_of_a_omega(l, m, a_array[i], abs(kerr_geo_orbital_frequency_circ(a_array[i], r_array[i]))) for i in range(a_array.shape[0])]
        else:
            amp_list = [[self.harmonicscpp.phase_of_a_omega(l, m, aa, abs(kerr_geo_orbital_frequency_circ(aa, rr))) for rr in r_array] for aa in a_array]
        return np.array(amp_list)
    
    # def __call__(self, int l, int m, double a, double r):
    #     return self.amplitude(l, m, a, r)*np.exp(1.j*self.phase(l, m, a, r))

    def key_check(self, int l, int m):
        return self.harmonicscpp.key_check(pair[int, int](l, m))
