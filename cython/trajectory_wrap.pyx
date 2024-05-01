from libcpp.vector cimport vector
import numpy as np
cimport numpy as np
from libc.string cimport memcpy
from libcpp.string cimport string
from libc.math cimport sin, cos, acos, abs, M_PI
from cython.operator import dereference
cimport openmp
import os

cdef unicode path_to_file = os.path.dirname(os.path.abspath(__file__))
cdef unicode default_trajectory_file = path_to_file + '/bhpwave/data/trajectory.txt'

cdef extern from "trajectory.hpp":
    cdef cppclass TrajectorySpline2D:
        TrajectorySpline2D(string filename) except +

        double time(double chi, double alpha)
        double phase(double chi, double alpha)
        double flux(double chi, double alpha)
        double flux_norm(double chi, double alpha)
        double orbital_alpha(double chi, double t)
        double orbital_alpha_derivative(double chi, double t)
        double orbital_frequency_time_derivative(double chi, double alpha)

        double time_of_a_omega(double a, double omega)
        double time_of_a_omega_derivative(double a, double omega)
        double phase_of_a_omega(double a, double omega)
        double phase_of_a_omega_derivative(double a, double omega)
        double flux_of_a_omega(double a, double omega)
        double orbital_frequency(double a, double t)
        double orbital_frequency_derivative(double a, double t)
        double orbital_frequency_time_derivative_of_a_omega(double a, double omega)
        double phase_of_a_time(double a, double t)

        double orbital_frequency_isco(double chi)
        double orbital_frequency_isco_of_a(double a)
        double min_orbital_frequency(double a)
        double max_orbital_frequency(double a)
        double max_orbital_radius(double a)
        double max_time_before_merger(double a)

        void flux_of_a_omega(double flux[], const double a[], const double omega[], int n, int num_threads)

    cdef cppclass InspiralContainer:
        InspiralContainer(int inspiralSteps)
        void setInspiralInitialConditions(double a, double massratio, double r0, double dt)
        void setTimeStep(int i, double alpha, double phase)
        void setTimeSteps(double* alpha, double* phase)

        const vector[double]& getAlpha() const
        const vector[double]& getPhase() const
        vector[double]& getAlphaNonConstRef()
        vector[double]& getPhaseNonConstRef()

        double getAlpha(int i)
        double getPhase(int i)

        double getTime(int i)
        double getFrequency(int i)
        double getRadius(int i)

        double getSpin()
        double getMassRatio()
        double getInitialRadius()
        double getInitialFrequency()
        int getSize()

    cdef cppclass InspiralGenerator:
        InspiralGenerator(TrajectorySpline2D &traj, int num_threads)
        # InspiralContainer computeInspiral(double a, double massratio, double r0, double dt, double T, int num_threads)
        void computeInitialConditions(double &chi, double &omega_i, double &alpha_i, double &t_i, double a, double massratio, double r0, double &T)
        int computeTimeStepNumber(double dt, double T)
        void computeInspiral(InspiralContainer &inspiral, double chi, double omega_i, double alpha_i, double t_i, double massratio, double dt, int num_threads)

######################################
# Define Useful Trajectory Functions #
######################################

OMEGA_MIN = 2.e-3
A_MAX = 0.9999

def kerr_geo_radius_circ(a, omega):
    return (abs(omega)*(1. - a*omega)/(omega**2))**(2./3.)

def kerr_geo_orbital_frequency_circ(a, r):
    v = 1./np.sqrt(r)
    return pow(v, 3)/(1 + a*pow(v, 3))

def kerr_isco_radius(a):
    sgnX = np.sign(a)
    z1 = 1 + pow(1 - a*a, 1./3.)*(pow(1 - a, 1./3.) + pow(1 + a, 1./3.))
    z2 = np.sqrt(3*a*a + z1*z1)

    return 3 + z2 - sgnX*np.sqrt((3. - z1)*(3. + z1 + 2.*z2))

def kerr_isco_frequency(a):
    rISCO = kerr_isco_radius(a)
    return kerr_geo_orbital_frequency_circ(a, rISCO)

def alpha_of_a_omega(a, omega):
    oISCO = kerr_isco_frequency(a)
    return alpha_of_omega_ISCO(omega, oISCO)

def alpha_of_omega_ISCO(omega, oISCO):
    return (abs(oISCO**(1./3.) - omega**(1./3.))/(oISCO**(1./3.) - OMEGA_MIN**(1./3.)))**(0.5)

def omega_of_a_alpha(a, alpha):
    oISCO = kerr_isco_frequency(a)
    return omega_of_alpha_ISCO(alpha, oISCO)

def omega_of_alpha_ISCO(alpha, oISCO):
    return pow(pow(oISCO, 1./3.) - pow(alpha, 2.)*(pow(oISCO, 1./3.) - pow(OMEGA_MIN, 1./3.)), 3.)

def chi_of_spin_subfunc(a):
    return pow(1. - a, 1./3.)

def chi_of_spin(a):
    return pow((chi_of_spin_subfunc(a) - chi_of_spin_subfunc(A_MAX))/(chi_of_spin_subfunc(-A_MAX) - chi_of_spin_subfunc(A_MAX)), 0.5)

def a_omega_to_chi_alpha(a, omega):
    chi = chi_of_spin(a)
    alpha = alpha_of_a_omega(a, omega)
    return (chi, alpha)

#####################################
# Define New Python Wrapped Classes #
#####################################

cdef class InspiralContainerWrapper:
    cdef InspiralContainer *inspiralcpp

    def __cinit__(self, int timeSteps=1):
        self.inspiralcpp = new InspiralContainer(timeSteps)

    def __dealloc__(self):
        del self.inspiralcpp

    def set_initial_conditions(self, double a, double epsilon, double r0, double dt):
        self.inspiralcpp.setInspiralInitialConditions(a, epsilon, r0, dt)

    def set_time_step(self, int i, double alpha, double phase):
        self.inspiralcpp.setTimeStep(i, alpha, phase)

    def set_time_steps(self, np.ndarray[ndim = 1, dtype = np.float64_t, mode="c"] alpha, np.ndarray[ndim = 1, dtype = np.float64_t, mode="c"] phase):
        self.inspiralcpp.setTimeSteps(&alpha[0], &phase[0])
    
    @property
    def timesteps(self):
        return self.inspiralcpp.getPhase().size()

    @property
    def size(self):
        return self.inspiralcpp.getPhase().size()

    @property
    def phase(self):
        cdef double[::1] arr = <double [:self.size]>self.inspiralcpp.getPhaseNonConstRef().data()
        return np.asarray(arr)

    @property
    def alpha(self):
        cdef double[::1] arr = <double [:self.inspiralcpp.getAlpha().size()]>self.inspiralcpp.getAlphaNonConstRef().data()
        return np.asarray(arr)

    @property
    def spin(self):
        return self.inspiralcpp.getSpin()
    
    @property
    def a(self):
        return self.inspiralcpp.getSpin()
    
    @property
    def massratio(self):
        return self.inspiralcpp.getMassRatio()
    
    @property
    def initialradius(self):
        return self.inspiralcpp.getInitialRadius()

    @property
    def initialfrequency(self):
        return self.inspiralcpp.getInitialFrequency()

    @property
    def iscofrequency(self):
        return kerr_isco_frequency(self.a)

    @property
    def dt(self):
        return self.inspiralcpp.getTime(1)

    @property
    def time(self):
        return np.arange(0, self.inspiralcpp.getSize())*self.dt

    @property
    def frequency(self):
        return omega_of_alpha_ISCO(self.alpha, self.iscofrequency)

    @property
    def radius(self):
        return kerr_geo_radius_circ(self.a, self.frequency)

cdef class InspiralGeneratorPy:
    cdef InspiralGenerator *inspiralcpp

    def __cinit__(self, TrajectoryDataPy traj, int num_threads=0):
        self.inspiralcpp = new InspiralGenerator(dereference(traj.trajcpp), num_threads)

    def __dealloc__(self):
        del self.inspiralcpp
        
    def __call__(self, double massratio, double a, double r0, double dt, double T, int num_threads=0):
        cdef double chi = 0.
        cdef double omega_i = 0.
        cdef double alpha_i = 0.
        cdef double t_i = 0.
        self.inspiralcpp.computeInitialConditions(chi, omega_i, alpha_i, t_i, a, massratio, r0, T)
        cdef int step_number = self.inspiralcpp.computeTimeStepNumber(dt, T)
        cdef InspiralContainerWrapper inspiral = InspiralContainerWrapper(step_number)
        inspiral.set_initial_conditions(a, massratio, r0, dt)
        self.inspiralcpp.computeInspiral(dereference(inspiral.inspiralcpp), chi, omega_i, alpha_i, t_i, massratio, dt, num_threads)
        # current design just copies since inspiral may be deleted after the function call. 
        # Not really efficient, but something more efficient requires initializing InspiralWrapper with pointers
        # to the memory types or numpy arrays that will live on
        # t = inspiral.time.copy()
        # r = inspiral.radius.copy()
        # phase = inspiral.phase.copy()
        # return (t, r, phase)
        return inspiral

import warnings

cdef class TrajectoryDataPy:
    cdef TrajectorySpline2D *trajcpp
    cdef bint dealloc_flag

    def __cinit__(self, unicode filename = default_trajectory_file, bint dealloc_flag = True):
        self.trajcpp = new TrajectorySpline2D(filename.encode())
        self.dealloc_flag = dealloc_flag

    def __dealloc__(self):
        if self.dealloc_flag:
            warnings.warn("Deallocating TrajectoryDataPy object", UserWarning)
        del self.trajcpp

    def time_to_merger(self, double a, double omega):
        return -self.trajcpp.time_of_a_omega(a, omega)

    def phase_to_merger(self, double a, double omega):
        return -self.trajcpp.phase_of_a_omega(a, omega)

    def phase_of_chi_alpha(self, double chi, double alpha):
        return self.trajcpp.phase(chi, alpha)

    def time_of_chi_alpha(self, double chi, double alpha):
        return self.trajcpp.time(chi, alpha)

    cdef flux_parallel(self, double a, np.ndarray[ndim = 1, dtype=np.float64_t, mode='c'] omega):
        cdef np.ndarray[ndim = 1, dtype=np.float64_t, mode='c'] flux = np.empty(omega.shape[0], dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype=np.float64_t, mode='c'] anp = a*np.ones(omega.shape[0], dtype=np.float64)
        self.trajcpp.flux_of_a_omega(&flux[0], &anp[0], &omega[0], omega.shape[0], 0)
        return flux

    def flux(self, double a, omega):
        if isinstance(omega, np.ndarray):
            return self.flux_parallel(a, omega)
        elif type(omega) is float or type(omega) is np.float64:
            return self.trajcpp.flux_of_a_omega(a, omega)
        else:
            raise TypeError("Frequency must be a float of numpy array")

    def orbital_frequency(self, double a, double t):
        return self.trajcpp.orbital_frequency(a, t)

    def orbital_alpha(self, double a, double t):
        return self.trajcpp.orbital_alpha(chi_of_spin(a), t)
    
    def phase(self, double a, double t):
        return -self.trajcpp.phase_of_a_time(a, t)

    def isco_frequency(self, double a):
        return self.trajcpp.orbital_frequency_isco_of_a(a)

    def min_orbital_frequency(self, double a):
        return self.trajcpp.min_orbital_frequency(a)

    def max_orbital_frequency(self, double a):
        return self.trajcpp.max_orbital_frequency(a)

    def max_time_before_merger(self, double a):
        return self.trajcpp.max_time_before_merger(a)
    
