from libcpp.vector cimport vector
from libc.string cimport memcpy 
import numpy as np
cimport numpy as np
from cython.operator import dereference

cdef extern from "swsh.hpp":
    vector[double] spin_weighted_spherical_harmonic(int, int, int, vector[double])
    double spin_weighted_spherical_harmonic(int, int, int, double)
    void spin_weighted_spherical_harmonic(double*, int, int, int, int, double*)

cdef spin_weighted_spherical_harmonic_vec_2(int s, int l, int m, np.ndarray[ndim=1, dtype=np.float64_t] theta):
    cdef np.ndarray[ndim=1, dtype=np.float64_t, mode='c'] vecpy = np.empty(theta.shape[0], dtype=np.float64)
    spin_weighted_spherical_harmonic(&vecpy[0], vecpy.shape[0], s, l, m, &theta[0])

    return vecpy

cdef spin_weighted_spherical_harmonic_double(int s, int l, int m, double theta):
    return spin_weighted_spherical_harmonic(s, l, m, theta)

def Yslm(int s, int l, int m, theta):
    if isinstance(theta, np.ndarray):
        return spin_weighted_spherical_harmonic_vec_2(s, l, m, theta)
    elif isinstance(theta, list):
        return spin_weighted_spherical_harmonic_vec_2(s, l, m, np.array(theta))
    else:
        return spin_weighted_spherical_harmonic_double(s, l, m, theta)