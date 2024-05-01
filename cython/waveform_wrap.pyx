from libcpp.vector cimport vector
import numpy as np
cimport numpy as np
from libcpp.string cimport string
from libcpp.complex cimport complex as cpp_complex

include "harmonic_wrap.pyx"
include "spline_wrap.pyx"

cdef extern from "waveform.hpp":
    double years_to_seconds(double years)
    double seconds_to_years(double seconds)
    double solar_mass_to_seconds(double mass)
    double scale_strain_amplitude(double mu, double dist)
    void sourceAngles(double &theta, double &phi, double qS, double phiS, double qK, double phiK);
    cpp_complex[double] polarization(double qS, double phiS, double qK, double phiK)

    cdef cppclass WaveformContainer:
        WaveformContainer(int timeSteps) except +
        WaveformContainer(double* plus, double *cross, int timeSteps) except +
        void setTimeStep(int i, double plus, double cross)
        void addTimeStep(int i, double plus, double cross)

        double* getPlusPointer()
        double* getCrossPointer()

        double getPlus(int i)
        double getCross(int i)
        int getSize()

    cdef cppclass WaveformHarmonicsContainer:
        WaveformHarmonicsContainer(int modeNum, int timeSteps) except +
        WaveformHarmonicsContainer(double *plus_ptr, double *cross_ptr, int modeNum, int timeSteps) except +
        void setTimeStep(int i, int j, double plus, double cross)
        void addTimeStep(int i, int j, double plus, double cross)
        void multiplyTimeStep(int i, int j, double plus, double cross)

        double* getPlusPointer()
        double* getCrossPointer()

        double getPlus(int i, int j)
        double getCross(int i, int j)

        int getSize()
        int getTimeSize()
        int getModeSize()

    cdef cppclass WaveformHarmonicOptions:
        WaveformHarmonicOptions()
        WaveformHarmonicOptions(double rescale, int num, int pad_output, int include_negative_m)

        cpp_complex[double] rescale
        int num_threads
        int pad_output
        int include_negative_m

    cdef cppclass WaveformHarmonicGenerator:
        WaveformHarmonicGenerator(HarmonicAmplitudes &Alm, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts)

        int computeTimeStepNumber(double dt, double T)
        WaveformContainer computeWaveformHarmonic(int l, int m, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts) except +
        void computeWaveformHarmonic(WaveformContainer &h, int l, int m, InspiralContainer &inspiral, double theta, double phi) except +
        void computeWaveformHarmonic(WaveformContainer &h, int l, int m, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts) except +

        WaveformContainer computeWaveformHarmonics(int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts)
        void computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts)
        void computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi, HarmonicOptions hOpts)
        void computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts) except +
        
        void computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts)
        void computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts) except +

        void computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi)
        void computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi)
        void computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, double theta, double phi)

        void computeWaveformHarmonics(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts)
        void computeWaveformHarmonics(WaveformHarmonicsContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts) except +
        
        void computeWaveformHarmonicsPhaseAmplitude(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts)
        void computeWaveformHarmonicsPhaseAmplitude(WaveformHarmonicsContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts) except +

        HarmonicSelector& getModeSelector()
        HarmonicModeContainer selectModes(InspiralContainer &inspiral, double theta)
        HarmonicModeContainer selectModes(InspiralContainer &inspiral, double theta, HarmonicOptions opts) except +
                
        WaveformHarmonicOptions getWaveformHarmonicOptions()
        HarmonicOptions getHarmonicOptions()

    cdef cppclass WaveformGenerator:
        WaveformGenerator(TrajectorySpline2D &traj, HarmonicAmplitudes &harm, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts)

        double convertTime(double t, double M)
        int computeTimeStepNumber(double dt, double T)
        int computeTimeStepNumber(double M, double mu, double a, double r0, double dt, double T)

        void computeWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T)
        void computeWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, WaveformHarmonicOptions wOpts)
        void computeWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions hOpts)
        void computeWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts)
        void computeWaveform(WaveformContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts)

        HarmonicModeContainer selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T)
        HarmonicModeContainer selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions opts)

        void computeWaveformSourceFrame(WaveformContainer &h, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double dt, double T)
        void computeWaveformSourceFrame(WaveformContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double dt, double T)

        void computeWaveform(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts)
        void computeWaveformSourceFrame(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double dt, double T)

        void computeWaveformPhaseAmplitude(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts)

        WaveformHarmonicOptions getWaveformHarmonicOptions()
        HarmonicOptions getHarmonicOptions()

cdef extern from "fourier.hpp":
    cdef cppclass WaveformFourierHarmonicGenerator:
        WaveformFourierHarmonicGenerator(HarmonicAmplitudes &Alm, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts)

        void computeWaveformFourierHarmonics(WaveformContainer &h, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, HarmonicOptions hOpts, int num_threads, double freq[], int fsamples)
        void computeWaveformFourierHarmonics(WaveformContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, int num_threads, double freq[], int fsamples)
        void computeWaveformFourierHarmonics(WaveformContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, int num_threads, double freq[], int fsamples)

        HarmonicSelector& getModeSelector()
        HarmonicModeContainer selectModes(InspiralContainer &inspiral, double theta)
        HarmonicModeContainer selectModes(InspiralContainer &inspiral, double theta, HarmonicOptions opts)

        WaveformHarmonicOptions getWaveformHarmonicOptions()
        HarmonicOptions getHarmonicOptions()

    cdef cppclass WaveformFourierGenerator:
        WaveformFourierGenerator(TrajectorySpline2D &traj, HarmonicAmplitudes &harm, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts)

        double convertTime(double t, double M)
        double convertFrequency(double f, double M)
        int computeFrequencyStepNumber(double df, double T)
        int computeTimeStepNumber(double dt, double T)
        int computeTimeStepNumber(double M, double mu, double a, double r0, double dt, double T)
        
        void computeFourierWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double freq[], double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts)
        void computeFourierWaveform(WaveformContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double freq[], double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts)
        
        void computeFourierWaveform(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double frequencies[], double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts)

        void computeFourierWaveformPhaseAmplitude(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double frequencies[], double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts)

        void computeFourierWaveformSourceFrame(WaveformContainer &h, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double freq[], double T)
        void computeFourierWaveformSourceFrame(WaveformContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double freq[], double T)

        void computeFourierWaveformSourceFrame(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double frequencies[], double T)

        HarmonicModeContainer selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double T)
        HarmonicModeContainer selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double T, HarmonicOptions opts)

        WaveformHarmonicOptions getWaveformHarmonicOptions()
        HarmonicOptions getHarmonicOptions()

# https://stackoverflow.com/questions/49400500/passing-1-or-2-d-numpy-array-to-c-throw-cython
cdef double* get_array_pointer(arr) except NULL:
    assert(arr.flags.c_contiguous) # if this isn't true, ravel will make a copy
    cdef double[::1] mview = arr.ravel()
    return &mview[0]

cdef class WaveformContainerNumpyWrapper:
    cdef WaveformContainer *hcpp

    def __cinit__(self, np.ndarray[ndim=1, dtype=np.float64_t, mode='c'] plus, np.ndarray[ndim=1, dtype=np.float64_t, mode='c'] cross):
        cdef steps = plus.shape[0]
        self.hcpp = new WaveformContainer(&plus[0], &cross[0], steps)

    def __dealloc__(self):
        del self.hcpp

    def plus_at_time_step(self, int i):
        return self.hcpp.getPlus(i)
    
    def cross_at_time_step(self, int i):
        return self.hcpp.getCross(i)

    @property
    def size(self):
        return self.hcpp.getSize()

    @property
    def plus(self):
        cdef double[::1] arr = <double [:self.size]>self.hcpp.getPlusPointer()
        return np.asarray(arr)

    @property
    def cross(self):
        cdef double[::1] arr = <double [:self.size]>self.hcpp.getCrossPointer()
        return np.asarray(arr)

cdef class WaveformHarmonicsContainerNumpyWrapper:
    cdef WaveformHarmonicsContainer *hcpp

    def __cinit__(self, np.ndarray[ndim=2, dtype=np.float64_t, mode='c'] plus, np.ndarray[ndim=2, dtype=np.float64_t, mode='c'] cross):
        cdef modeNum = plus.shape[0]
        cdef timeSteps = plus.shape[1]
        self.hcpp = new WaveformHarmonicsContainer(get_array_pointer(plus), get_array_pointer(cross), modeNum, timeSteps)

    def __dealloc__(self):
        del self.hcpp

    def plus_step(self, int i, int j):
        return self.hcpp.getPlus(i, j)
    
    def cross_step(self, int i, int j):
        return self.hcpp.getCross(i, j)

    @property
    def size(self):
        return (self.hcpp.getModeSize(), self.hcpp.getTimeSize())

    @property
    def plus(self):
        cdef double[:,::1] arr = <double [:self.size[0], :self.size[1]]>self.hcpp.getPlusPointer()
        return np.asarray(arr)

    @property
    def cross(self):
        cdef double[:,::1] arr = <double [:self.size[0], :self.size[1]]>self.hcpp.getCrossPointer()
        return np.asarray(arr)

cdef class WaveformContainerWrapper:
    cdef WaveformContainer *hcpp

    def __cinit__(self, int time_steps):
        self.hcpp = new WaveformContainer(time_steps)

    def __dealloc__(self):
        del self.hcpp

    def plus_at_time_step(self, int i):
        return self.hcpp.getPlus(i)
    
    def cross_at_time_step(self, int i):
        return self.hcpp.getPlus(i)

    @property
    def size(self):
        return self.hcpp.getSize()

    @property
    def plus(self):
        cdef double[::1] arr = <double [:self.size]>self.hcpp.getPlusPointer()
        return np.asarray(arr)

    @property
    def cross(self):
        cdef double[::1] arr = <double [:self.size]>self.hcpp.getCrossPointer()
        return np.asarray(arr)

cdef class WaveformGeneratorPy:
    cdef WaveformGenerator *hcpp

    def __cinit__(self, TrajectoryDataPy traj, HarmonicAmplitudesPy Alm, dict harmonic_kwargs = {}, dict waveform_kwargs = {}):
        cdef WaveformHarmonicOptions wOpts
        cdef HarmonicOptions hOpts
        
        if "eps" in harmonic_kwargs.keys():
            hOpts.epsilon = harmonic_kwargs["eps"]
        if "max_samples" in harmonic_kwargs.keys():
            hOpts.max_samples = harmonic_kwargs["max_samples"]
        
        if "num_threads" in waveform_kwargs.keys():
            wOpts.num_threads = waveform_kwargs["num_threads"]
        if "pad_output" in waveform_kwargs.keys():
            wOpts.pad_output = waveform_kwargs["pad_output"]
        if "include_negative_m" in waveform_kwargs.keys():
            wOpts.include_negative_m = waveform_kwargs["include_negative_m"]

        self.hcpp = new WaveformGenerator(dereference(traj.trajcpp), dereference(Alm.harmonicscpp), hOpts, wOpts)

    def __dealloc__(self):
        del self.hcpp

    def time_step_number(self, double M, double mu, double a, double r0, double dt, double T):
        return self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)

    def select_modes(self, double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, pad_nmodes = False, *args, **kwargs):
        cdef HarmonicOptions hOpts
        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        cdef HarmonicModeContainer modescpp = self.hcpp.selectModes(M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, dt, T, hOpts)
        cdef HarmonicModeContainerWrapper modeWrap = HarmonicModeContainerWrapper()
        modeWrap.wrap(modescpp)

        if pad_nmodes:
            select_modes = list(zip(modeWrap.lmodes, modeWrap.mmodes, 0*modeWrap.mmodes))
        else:
            select_modes = list(zip(modeWrap.lmodes, modeWrap.mmodes))

        return select_modes

    def waveform_harmonics(self, int[::1] l, int[::1] m, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(timeSteps, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(timeSteps, dtype=np.float64)
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross)
        
        self.hcpp.computeWaveform(dereference(h.hcpp), &l[0], &m[0], l.shape[0], M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, hOpts, wOpts)
        if return_list:
            return [plus, cross]

        # cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] waveform = np.empty(timeSteps, dtype=np.complex128)
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] waveform = -1.j*cross
        # waveform = -1.j*cross
        waveform += plus

        return waveform

    def waveform_harmonics_grid(self, int[::1] l, int[::1] m, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int timeSteps
        cdef int modeNum = l.shape[0]
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] plus = np.zeros((modeNum, timeSteps), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] cross = np.zeros((modeNum, timeSteps), dtype=np.float64)
        cdef WaveformHarmonicsContainerNumpyWrapper h = WaveformHarmonicsContainerNumpyWrapper(plus, cross)
        
        self.hcpp.computeWaveform(dereference(h.hcpp), &l[0], &m[0], modeNum, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, hOpts, wOpts)
        if return_list:
            return [plus, cross]

        # cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] waveform = np.empty((modeNum, timeSteps), dtype=np.complex128)
        cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] waveform = -1.j*cross
        waveform += plus

        return waveform

    def waveform_select_harmonics_grid(self, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef HarmonicModeContainer modescpp = self.hcpp.selectModes(M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, dt, T, hOpts)
        cdef HarmonicModeContainerWrapper modeWrap = HarmonicModeContainerWrapper()
        modeWrap.wrap(modescpp)
        # cdef np.ndarray[ndim = 1, dtype = np.int32_t, mode='c'] l = modeWrap.lmodes
        # cdef np.ndarray[ndim = 1, dtype = np.int32_t, mode='c'] m = modeWrap.mmodes
        cdef int[::1] l = modeWrap.lmodes.data
        cdef int[::1] m = modeWrap.mmodes.data
        cdef int modeNum = len(l)

        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] plus = np.zeros((modeNum, timeSteps), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] cross = np.zeros((modeNum, timeSteps), dtype=np.float64)
        cdef WaveformHarmonicsContainerNumpyWrapper h = WaveformHarmonicsContainerNumpyWrapper(plus, cross)
        
        self.hcpp.computeWaveform(dereference(h.hcpp), &l[0], &m[0], modeNum, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, hOpts, wOpts)
        if return_list:
            return [plus, cross]

        # cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] waveform = np.empty((modeNum, timeSteps), dtype=np.complex128)
        cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] waveform = -1.j*cross
        waveform += plus

        return waveform

    def waveform_harmonics_phase_amplitude(self, int[::1] l, int[::1] m, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, *args, **kwargs):
        cdef int timeSteps
        cdef int modeNum = l.shape[0]
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] amp = np.zeros((modeNum, timeSteps), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] phase = np.zeros((modeNum, timeSteps), dtype=np.float64)
        cdef WaveformHarmonicsContainerNumpyWrapper h = WaveformHarmonicsContainerNumpyWrapper(amp, phase)
        
        self.hcpp.computeWaveformPhaseAmplitude(dereference(h.hcpp), &l[0], &m[0], modeNum, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, hOpts, wOpts)
        return [amp, phase]

    def waveform_select_harmonics_phase_amplitude(self, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, *args, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef HarmonicModeContainer modescpp = self.hcpp.selectModes(M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, dt, T, hOpts)
        cdef HarmonicModeContainerWrapper modeWrap = HarmonicModeContainerWrapper()
        modeWrap.wrap(modescpp)

        cdef int[::1] l = modeWrap.lmodes.data
        cdef int[::1] m = modeWrap.mmodes.data
        cdef int modeNum = len(l)

        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] amp = np.zeros((2*modeNum, timeSteps), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] phase = np.zeros((2*modeNum, timeSteps), dtype=np.float64)
        cdef WaveformHarmonicsContainerNumpyWrapper h = WaveformHarmonicsContainerNumpyWrapper(amp, phase)
        
        self.hcpp.computeWaveformPhaseAmplitude(dereference(h.hcpp), &l[0], &m[0], modeNum, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, hOpts, wOpts)
        return [amp, phase]
    
    def waveform(self, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list = False, *args, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]
        # cdef WaveformContainerWrapper h = WaveformContainerWrapper(timeSteps)

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(timeSteps, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(timeSteps, dtype=np.float64)
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross)

        self.hcpp.computeWaveform(dereference(h.hcpp), M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, hOpts, wOpts)

        if return_list:
            return [plus, cross]
    
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] waveform = np.empty(timeSteps, dtype=np.complex128)
        waveform = -1.j*cross
        waveform += plus

        return waveform

    def waveform_harmonics_source_frame(self, int[::1] l, int[::1] m, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(timeSteps, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(timeSteps, dtype=np.float64)
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross)
        
        self.hcpp.computeWaveformSourceFrame(dereference(h.hcpp), &l[0], &m[0], l.shape[0], M, mu, a, r0, theta, phi, Phi_phi0, dt, T)
        if return_list:
            return [plus, cross]

        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] waveform = np.empty(timeSteps, dtype=np.complex128)
        waveform = -1.j*cross
        waveform += plus

        return waveform

    def waveform_source_frame(self, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(timeSteps, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(timeSteps, dtype=np.float64)
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross)

        self.hcpp.computeWaveformSourceFrame(dereference(h.hcpp), M, mu, a, r0, theta, phi, Phi_phi0, dt, T)
        if return_list:
            return [plus, cross]

        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] waveform = np.empty(timeSteps, dtype=np.complex128)
        waveform = -1.j*cross
        waveform += plus

        return waveform

cdef class WaveformFourierGeneratorPy:
    cdef WaveformFourierGenerator *hcpp

    def __cinit__(self, TrajectoryDataPy traj, HarmonicAmplitudesPy Alm, dict harmonic_kwargs = {}, dict waveform_kwargs = {}):
        cdef WaveformHarmonicOptions wOpts
        cdef HarmonicOptions hOpts
        
        if "eps" in harmonic_kwargs.keys():
            hOpts.epsilon = harmonic_kwargs["eps"]
        if "max_samples" in harmonic_kwargs.keys():
            hOpts.max_samples = harmonic_kwargs["max_samples"]
        
        if "num_threads" in waveform_kwargs.keys():
            wOpts.num_threads = waveform_kwargs["num_threads"]
        if "pad_output" in waveform_kwargs.keys():
            wOpts.pad_output = waveform_kwargs["pad_output"]
        if "include_negative_m" in waveform_kwargs.keys():
            wOpts.include_negative_m = waveform_kwargs["include_negative_m"]

        self.hcpp = new WaveformFourierGenerator(dereference(traj.trajcpp), dereference(Alm.harmonicscpp), hOpts, wOpts)

    def __dealloc__(self):
        del self.hcpp

    def step_number(self, double dt, double T):
        return self.hcpp.computeFrequencyStepNumber(dt, T)

    def select_modes(self, double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double T, pad_nmodes = False, *args, **kwargs):
        cdef HarmonicOptions hOpts
        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        cdef HarmonicModeContainer modescpp = self.hcpp.selectModes(M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, T, hOpts)
        cdef HarmonicModeContainerWrapper modeWrap = HarmonicModeContainerWrapper()
        modeWrap.wrap(modescpp)

        if pad_nmodes:
            select_modes = list(zip(modeWrap.lmodes, modeWrap.mmodes, 0*modeWrap.mmodes))
        else:
            select_modes = list(zip(modeWrap.lmodes, modeWrap.mmodes))

        return select_modes

    def waveform_harmonics(self, int[::1] l, int[::1] m, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double[::1] frequencies, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int steps
        cdef int steps_no_zero
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        # if the user specifies a value for dt, then pick frequencies to align with DFT of TD signal
        if dt > 0.:
            if pad_output:
                timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
            else:
                timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
            frequencies = np.fft.rfftfreq(timeSteps, d=dt)

        steps = len(frequencies)
        if frequencies[0] == 0.:
            frequencies = frequencies[1:] # throw away zero frequencies
        steps_no_zero = len(frequencies) # number of steps without zero frequency

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(2*steps_no_zero, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(2*steps_no_zero, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] plusComplex = np.zeros(steps, dtype=np.complex128)
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] crossComplex = np.zeros(steps, dtype=np.complex128)

        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross)

        self.hcpp.computeFourierWaveform(dereference(h.hcpp), &l[0], &m[0], l.shape[0], M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, &frequencies[0], T, hOpts, wOpts)
        plusComplex[steps-steps_no_zero:] = 1.j*np.flip(plus[-steps_no_zero:])
        plusComplex[steps-steps_no_zero:] += plus[:steps_no_zero]
        crossComplex[steps-steps_no_zero:] = 1.j*np.flip(cross[-steps_no_zero:])
        crossComplex[steps-steps_no_zero:] += cross[:steps_no_zero]
        
        return [plusComplex, crossComplex]

    def waveform_harmonics_grid(self, int[::1] l, int[::1] m, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double[::1] frequencies, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int steps
        cdef int steps_no_zero
        cdef int timeSteps
        cdef int modeNum
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        # if the user specifies a value for dt, then pick frequencies to align with DFT of TD signal
        if dt > 0.:
            if pad_output:
                timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
            else:
                timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
            frequencies = np.fft.rfftfreq(timeSteps, d=dt)

        steps = len(frequencies)
        modeNum = len(l)
        if frequencies[0] == 0.:
            frequencies = frequencies[1:] # throw away zero frequencies
        steps_no_zero = len(frequencies) # number of steps without zero frequency

        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] plus = np.zeros((modeNum, 2*steps_no_zero), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] cross = np.zeros((modeNum, 2*steps_no_zero), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] plusComplex = np.zeros((modeNum, steps), dtype=np.complex128)
        cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] crossComplex = np.zeros((modeNum, steps), dtype=np.complex128)
        cdef WaveformHarmonicsContainerNumpyWrapper h = WaveformHarmonicsContainerNumpyWrapper(plus, cross)

        self.hcpp.computeFourierWaveform(dereference(h.hcpp), &l[0], &m[0], l.shape[0], M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, &frequencies[0], T, hOpts, wOpts)
        plusComplex[:, steps-steps_no_zero:] = 1.j*np.flip(plus[:, -steps_no_zero:], axis = 1)
        plusComplex[:, steps-steps_no_zero:] += plus[:, :steps_no_zero]
        crossComplex[:, steps-steps_no_zero:] = 1.j*np.flip(cross[:, -steps_no_zero:], axis = 1)
        crossComplex[:, steps-steps_no_zero:] += cross[:, :steps_no_zero]
        
        return [plusComplex, crossComplex]

    def waveform_select_harmonics_grid(self, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double[::1] frequencies, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int steps
        cdef int steps_no_zero
        cdef int timeSteps
        cdef int modeNum
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        # if the user specifies a value for dt, then pick frequencies to align with DFT of TD signal
        if dt > 0.:
            if pad_output:
                timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
            else:
                timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
            frequencies = np.fft.rfftfreq(timeSteps, d=dt)

        cdef HarmonicModeContainer modescpp = self.hcpp.selectModes(M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, T, hOpts)
        cdef HarmonicModeContainerWrapper modeWrap = HarmonicModeContainerWrapper()
        modeWrap.wrap(modescpp)

        cdef int[::1] l = modeWrap.lmodes.data
        cdef int[::1] m = modeWrap.mmodes.data

        steps = len(frequencies)
        modeNum = len(l)
        if frequencies[0] == 0.:
            frequencies = frequencies[1:] # throw away zero frequencies
        steps_no_zero = len(frequencies) # number of steps without zero frequency

        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] plus = np.zeros((modeNum, 2*steps_no_zero), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] cross = np.zeros((modeNum, 2*steps_no_zero), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] plusComplex = np.zeros((modeNum, steps), dtype=np.complex128)
        cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] crossComplex = np.zeros((modeNum, steps), dtype=np.complex128)

        cdef WaveformHarmonicsContainerNumpyWrapper h = WaveformHarmonicsContainerNumpyWrapper(plus, cross)

        self.hcpp.computeFourierWaveform(dereference(h.hcpp), &l[0], &m[0], l.shape[0], M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, &frequencies[0], T, hOpts, wOpts)
        plusComplex[:, steps-steps_no_zero:] = 1.j*np.flip(plus[:, -steps_no_zero:], axis = 1)
        plusComplex[:, steps-steps_no_zero:] += plus[:, :steps_no_zero]
        crossComplex[:, steps-steps_no_zero:] = 1.j*np.flip(cross[:, -steps_no_zero:], axis = 1)
        crossComplex[:, steps-steps_no_zero:] += cross[:, :steps_no_zero]
        
        return [plusComplex, crossComplex]

    def waveform_harmonics_phase_amplitude(self, int[::1] l, int[::1] m, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double[::1] frequencies, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int steps
        cdef int steps_no_zero
        cdef int timeSteps
        cdef int modeNum
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        # if the user specifies a value for dt, then pick frequencies to align with DFT of TD signal
        if dt > 0.:
            if pad_output:
                timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
            else:
                timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
            frequencies = np.fft.rfftfreq(timeSteps, d=dt)

        steps = len(frequencies)
        modeNum = len(l)
        if frequencies[0] == 0.:
            frequencies = frequencies[1:] # throw away zero frequencies
        steps_no_zero = len(frequencies) # number of steps without zero frequency

        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] ampNoZero = np.zeros((2*modeNum, steps_no_zero), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] phaseNoZero = np.zeros((2*modeNum, steps_no_zero), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] amp = np.zeros((2*modeNum, steps), dtype=np.complex128)
        cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] phase = np.zeros((2*modeNum, steps), dtype=np.complex128)
        cdef WaveformHarmonicsContainerNumpyWrapper h = WaveformHarmonicsContainerNumpyWrapper(ampNoZero, phaseNoZero)

        self.hcpp.computeFourierWaveformPhaseAmplitude(dereference(h.hcpp), &l[0], &m[0], l.shape[0], M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, &frequencies[0], T, hOpts, wOpts)
        amp[:, steps-steps_no_zero:] = ampNoZero
        phase[:, steps-steps_no_zero:] = phaseNoZero
        
        return [amp, phase]

    def waveform_select_harmonics_phase_amplitude(self, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double[::1] frequencies, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int steps
        cdef int steps_no_zero
        cdef int timeSteps
        cdef int modeNum
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        # if the user specifies a value for dt, then pick frequencies to align with DFT of TD signal
        if dt > 0.:
            if pad_output:
                timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
            else:
                timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
            frequencies = np.fft.rfftfreq(timeSteps, d=dt)

        cdef HarmonicModeContainer modescpp = self.hcpp.selectModes(M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, T, hOpts)
        cdef HarmonicModeContainerWrapper modeWrap = HarmonicModeContainerWrapper()
        modeWrap.wrap(modescpp)

        cdef int[::1] l = modeWrap.lmodes.data
        cdef int[::1] m = modeWrap.mmodes.data

        steps = len(frequencies)
        modeNum = len(l)
        if frequencies[0] == 0.:
            frequencies = frequencies[1:] # throw away zero frequencies
        steps_no_zero = len(frequencies) # number of steps without zero frequency

        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] ampNoZero = np.zeros((2*modeNum, steps_no_zero), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] phaseNoZero = np.zeros((2*modeNum, steps_no_zero), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] amp = np.zeros((2*modeNum, steps), dtype=np.complex128)
        cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] phase = np.zeros((2*modeNum, steps), dtype=np.complex128)
        cdef WaveformHarmonicsContainerNumpyWrapper h = WaveformHarmonicsContainerNumpyWrapper(ampNoZero, phaseNoZero)

        self.hcpp.computeFourierWaveformPhaseAmplitude(dereference(h.hcpp), &l[0], &m[0], l.shape[0], M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, &frequencies[0], T, hOpts, wOpts)
        amp[:, steps-steps_no_zero:] = ampNoZero
        phase[:, steps-steps_no_zero:] = phaseNoZero
        
        return [amp, phase]

    def waveform(self, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double[::1] frequencies, double dt, double T, bint pad_output = False, bint return_list = False, *args, **kwargs):
        cdef int steps
        cdef int steps_no_zero
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        # if the user specifies a value for dt, then pick frequencies to align with DFT of TD signal
        if dt > 0.:
            if pad_output:
                timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
            else:
                timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
            frequencies = np.fft.rfftfreq(timeSteps, d=dt)

        steps = len(frequencies)
        if frequencies[0] == 0.:
            frequencies = frequencies[1:] # throw away zero frequencies
        steps_no_zero = len(frequencies) # number of steps without zero frequency

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(2*steps_no_zero, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(2*steps_no_zero, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] plusComplex = np.zeros(steps, dtype=np.complex128)
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] crossComplex = np.zeros(steps, dtype=np.complex128)
        
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross)

        self.hcpp.computeFourierWaveform(dereference(h.hcpp), M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, &frequencies[0], T, hOpts, wOpts)
        plusComplex[steps-steps_no_zero:] = 1.j*np.flip(plus[-steps_no_zero:])
        plusComplex[steps-steps_no_zero:] += plus[:steps_no_zero]
        crossComplex[steps-steps_no_zero:] = 1.j*np.flip(cross[-steps_no_zero:])
        crossComplex[steps-steps_no_zero:] += cross[:steps_no_zero]

        return [plusComplex, crossComplex]
    
    def waveform_harmonics_source_frame(self, int[::1] l, int[::1] m, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double[::1] frequencies, double dt, double T, bint pad_output = False, bint return_list = False, *args, **kwargs):
        cdef int steps
        cdef int steps_no_zero
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()

        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        # if the user specifies a value for dt, then pick frequencies to align with DFT of TD signal
        if dt > 0.:
            if pad_output:
                timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
            else:
                timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
            frequencies = np.fft.rfftfreq(timeSteps, d=dt)

        steps = len(frequencies)
        if frequencies[0] == 0.:
            frequencies = frequencies[1:] # throw away zero frequencies
        steps_no_zero = len(frequencies) # number of steps without zero frequency

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(2*steps_no_zero, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(2*steps_no_zero, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] plusComplex = np.zeros(steps, dtype=np.complex128)
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] crossComplex = np.zeros(steps, dtype=np.complex128)
        
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross)

        self.hcpp.computeFourierWaveformSourceFrame(dereference(h.hcpp), &l[0], &m[0], l.shape[0], M, mu, a, r0, theta, phi, Phi_phi0, &frequencies[0], T)
        plusComplex[steps-steps_no_zero:] = 1.j*np.flip(plus[-steps_no_zero:])
        plusComplex[steps-steps_no_zero:] += plus[:steps_no_zero]
        crossComplex[steps-steps_no_zero:] = 1.j*np.flip(cross[-steps_no_zero:])
        crossComplex[steps-steps_no_zero:] += cross[:steps_no_zero]

        return [plusComplex, crossComplex]
    
    def waveform_source_frame(self, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double[::1] frequencies, double dt, double T, bint pad_output = False, bint return_list = False, *args, **kwargs):
        cdef int steps
        cdef int steps_no_zero
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        
        if "pad_output" in kwargs.keys():
            pad_output = kwargs["pad_output"]
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        # if the user specifies a value for dt, then pick frequencies to align with DFT of TD signal
        if dt > 0.:
            if pad_output:
                timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
            else:
                timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
            frequencies = np.fft.rfftfreq(timeSteps, d=dt)

        steps = len(frequencies)
        if frequencies[0] == 0.:
            frequencies = frequencies[1:] # throw away zero frequencies
        steps_no_zero = len(frequencies) # number of steps without zero frequency

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(2*steps_no_zero, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(2*steps_no_zero, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] plusComplex = np.zeros(steps, dtype=np.complex128)
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] crossComplex = np.zeros(steps, dtype=np.complex128)
        
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross)

        self.hcpp.computeFourierWaveformSourceFrame(dereference(h.hcpp), M, mu, a, r0, theta, phi, Phi_phi0, &frequencies[0], T)
        plusComplex[steps-steps_no_zero:] = 1.j*np.flip(plus[-steps_no_zero:])
        plusComplex[steps-steps_no_zero:] += plus[:steps_no_zero]
        crossComplex[steps-steps_no_zero:] = 1.j*np.flip(cross[-steps_no_zero:])
        crossComplex[steps-steps_no_zero:] += cross[:steps_no_zero]

        return [plusComplex, crossComplex]

# Custom
cdef class WaveformHarmonicGeneratorPy:
    cdef WaveformHarmonicGenerator *hcpp

    def __cinit__(self, HarmonicAmplitudesPy Alm, dict harmonic_kwargs = {}, dict waveform_kwargs = {}):
        cdef WaveformHarmonicOptions wOpts
        cdef HarmonicOptions hOpts
        
        if "eps" in harmonic_kwargs.keys():
            hOpts.epsilon = harmonic_kwargs["eps"]
        if "max_samples" in harmonic_kwargs.keys():
            hOpts.max_samples = harmonic_kwargs["max_samples"]
        
        if "num_threads" in waveform_kwargs.keys():
            wOpts.num_threads = waveform_kwargs["num_threads"]
        if "pad_output" in waveform_kwargs.keys():
            wOpts.pad_output = waveform_kwargs["pad_output"]
        if "include_negative_m" in waveform_kwargs.keys():
            wOpts.include_negative_m = waveform_kwargs["include_negative_m"]

        self.hcpp = new WaveformHarmonicGenerator(dereference(Alm.harmonicscpp), hOpts, wOpts)

    def __dealloc__(self):
        del self.hcpp

    def select_modes(self, InspiralContainerWrapper inspiral, double qS, double phiS, double qK, double phiK, pad_nmodes = False, *args, **kwargs):
        cdef HarmonicOptions hOpts
        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        cdef double theta
        cdef double phi
        sourceAngles(theta, phi, qS, phiS, qK, phiK)

        cdef HarmonicModeContainer modescpp = self.hcpp.selectModes(dereference(inspiral.inspiralcpp), theta, hOpts)
        cdef HarmonicModeContainerWrapper modeWrap = HarmonicModeContainerWrapper()
        modeWrap.wrap(modescpp)

        if pad_nmodes:
            select_modes = list(zip(modeWrap.lmodes, modeWrap.mmodes, 0*modeWrap.mmodes))
        else:
            select_modes = list(zip(modeWrap.lmodes, modeWrap.mmodes))

        return select_modes

    def waveform_harmonics(self, int[::1] l, int[::1] m, double mu, InspiralContainerWrapper inspiral, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = inspiral.timesteps
        
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(timeSteps, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(timeSteps, dtype=np.float64)
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross)

        cdef double theta
        cdef double phi
        sourceAngles(theta, phi, qS, phiS, qK, phiK)
        wOpts.rescale = polarization(qS, phiS, qK, phiK)
        wOpts.rescale *= scale_strain_amplitude(mu, dist)
        
        self.hcpp.computeWaveformHarmonics(dereference(h.hcpp), &l[0], &m[0], l.shape[0], dereference(inspiral.inspiralcpp), theta, phi - Phi_phi0, wOpts)
        if return_list:
            return [plus, cross]

        # cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] waveform = np.empty(timeSteps, dtype=np.complex128)
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] waveform = -1.j*cross
        # waveform = -1.j*cross
        waveform += plus

        return waveform

    def waveform_harmonics_grid(self, int[::1] l, int[::1] m, double mu, InspiralContainerWrapper inspiral, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int timeSteps
        cdef int modeNum = l.shape[0]
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = inspiral.timesteps
        
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] plus = np.zeros((modeNum, timeSteps), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] cross = np.zeros((modeNum, timeSteps), dtype=np.float64)
        cdef WaveformHarmonicsContainerNumpyWrapper h = WaveformHarmonicsContainerNumpyWrapper(plus, cross)
        
        cdef double theta
        cdef double phi
        sourceAngles(theta, phi, qS, phiS, qK, phiK)
        wOpts.rescale = polarization(qS, phiS, qK, phiK)
        wOpts.rescale *= scale_strain_amplitude(mu, dist)

        self.hcpp.computeWaveformHarmonics(dereference(h.hcpp), &l[0], &m[0], modeNum, dereference(inspiral.inspiralcpp), theta, phi - Phi_phi0, wOpts)
        if return_list:
            return [plus, cross]

        cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] waveform = -1.j*cross
        waveform += plus

        return waveform

    def waveform_select_harmonics_grid(self, double mu, InspiralContainerWrapper inspiral, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = inspiral.timesteps
        
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef double theta
        cdef double phi
        sourceAngles(theta, phi, qS, phiS, qK, phiK)

        cdef HarmonicModeContainer modescpp = self.hcpp.selectModes(dereference(inspiral.inspiralcpp), theta, hOpts)
        cdef HarmonicModeContainerWrapper modeWrap = HarmonicModeContainerWrapper()
        modeWrap.wrap(modescpp)
        # cdef np.ndarray[ndim = 1, dtype = np.int32_t, mode='c'] l = modeWrap.lmodes
        # cdef np.ndarray[ndim = 1, dtype = np.int32_t, mode='c'] m = modeWrap.mmodes
        cdef int[::1] l = modeWrap.lmodes.data
        cdef int[::1] m = modeWrap.mmodes.data
        cdef int modeNum = len(l)

        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] plus = np.zeros((modeNum, timeSteps), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] cross = np.zeros((modeNum, timeSteps), dtype=np.float64)
        cdef WaveformHarmonicsContainerNumpyWrapper h = WaveformHarmonicsContainerNumpyWrapper(plus, cross)
        
        wOpts.rescale = polarization(qS, phiS, qK, phiK)
        wOpts.rescale *= scale_strain_amplitude(mu, dist)

        self.hcpp.computeWaveformHarmonics(dereference(h.hcpp), &l[0], &m[0], modeNum, dereference(inspiral.inspiralcpp), theta, phi - Phi_phi0, wOpts)
        if return_list:
            return [plus, cross]

        # cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] waveform = np.empty((modeNum, timeSteps), dtype=np.complex128)
        cdef np.ndarray[ndim = 2, dtype = np.complex128_t, mode='c'] waveform = -1.j*cross
        waveform += plus

        return waveform

    def waveform_harmonics_phase_amplitude(self, int[::1] l, int[::1] m, double mu, InspiralContainerWrapper inspiral, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, *args, **kwargs):
        cdef int timeSteps
        cdef int modeNum = l.shape[0]
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = inspiral.timesteps

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] amp = np.zeros((modeNum, timeSteps), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] phase = np.zeros((modeNum, timeSteps), dtype=np.float64)
        cdef WaveformHarmonicsContainerNumpyWrapper h = WaveformHarmonicsContainerNumpyWrapper(amp, phase)
        
        cdef double theta
        cdef double phi
        sourceAngles(theta, phi, qS, phiS, qK, phiK)
        wOpts.rescale = polarization(qS, phiS, qK, phiK)
        wOpts.rescale *= scale_strain_amplitude(mu, dist)

        self.hcpp.computeWaveformHarmonicsPhaseAmplitude(dereference(h.hcpp), &l[0], &m[0], modeNum, dereference(inspiral.inspiralcpp), theta, phi - Phi_phi0, wOpts)
        return [amp, phase]

    def waveform_select_harmonics_phase_amplitude(self, double mu, InspiralContainerWrapper inspiral, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, *args, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = inspiral.timesteps

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef double theta
        cdef double phi
        sourceAngles(theta, phi, qS, phiS, qK, phiK)

        cdef HarmonicModeContainer modescpp = self.hcpp.selectModes(dereference(inspiral.inspiralcpp), theta, hOpts)
        cdef HarmonicModeContainerWrapper modeWrap = HarmonicModeContainerWrapper()
        modeWrap.wrap(modescpp)

        cdef int[::1] l = modeWrap.lmodes.data
        cdef int[::1] m = modeWrap.mmodes.data
        cdef int modeNum = len(l)

        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] amp = np.zeros((2*modeNum, timeSteps), dtype=np.float64)
        cdef np.ndarray[ndim = 2, dtype = np.float64_t, mode='c'] phase = np.zeros((2*modeNum, timeSteps), dtype=np.float64)
        cdef WaveformHarmonicsContainerNumpyWrapper h = WaveformHarmonicsContainerNumpyWrapper(amp, phase)
        
        wOpts.rescale = polarization(qS, phiS, qK, phiK)
        wOpts.rescale *= scale_strain_amplitude(mu, dist)

        self.hcpp.computeWaveformHarmonicsPhaseAmplitude(dereference(h.hcpp), &l[0], &m[0], modeNum, dereference(inspiral.inspiralcpp), theta, phi - Phi_phi0, wOpts)
        return [amp, phase]
    
    def waveform(self, double mu, InspiralContainerWrapper inspiral, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list = False, *args, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = inspiral.timesteps

        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]
        # cdef WaveformContainerWrapper h = WaveformContainerWrapper(timeSteps)

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(timeSteps, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(timeSteps, dtype=np.float64)
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross)

        cdef double theta
        cdef double phi
        sourceAngles(theta, phi, qS, phiS, qK, phiK)

        cdef HarmonicModeContainer modescpp = self.hcpp.selectModes(dereference(inspiral.inspiralcpp), theta, hOpts)
        cdef HarmonicModeContainerWrapper modeWrap = HarmonicModeContainerWrapper()
        modeWrap.wrap(modescpp)

        cdef int[::1] l = modeWrap.lmodes.data
        cdef int[::1] m = modeWrap.mmodes.data
        cdef int modeNum = len(l)
        
        wOpts.rescale = polarization(qS, phiS, qK, phiK)
        wOpts.rescale *= scale_strain_amplitude(mu, dist)

        self.hcpp.computeWaveformHarmonics(dereference(h.hcpp), &l[0], &m[0], modeNum, dereference(inspiral.inspiralcpp), theta, phi - Phi_phi0, wOpts)

        if return_list:
            return [plus, cross]
    
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] waveform = np.empty(timeSteps, dtype=np.complex128)
        waveform = -1.j*cross
        waveform += plus

        return waveform

    def waveform_harmonics_source_frame(self, int[::1] l, int[::1] m, InspiralContainerWrapper inspiral, double theta, double phi, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = inspiral.timesteps
        
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(timeSteps, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(timeSteps, dtype=np.float64)
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross)

        self.hcpp.computeWaveformHarmonics(dereference(h.hcpp), &l[0], &m[0], l.shape[0], dereference(inspiral.inspiralcpp), theta, phi - Phi_phi0)
        if return_list:
            return [plus, cross]

        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] waveform = np.empty(timeSteps, dtype=np.complex128)
        waveform = -1.j*cross
        waveform += plus

        return waveform

    def waveform_source_frame(self, InspiralContainerWrapper inspiral, double theta, double phi, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list=False, *args, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = inspiral.timesteps
        
        if "return_list" in kwargs.keys():
            return_list = kwargs["return_list"]

        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        if "include_negative_m" in kwargs.keys():
            wOpts.include_negative_m = kwargs["include_negative_m"]

        cdef HarmonicModeContainer modescpp = self.hcpp.selectModes(dereference(inspiral.inspiralcpp), theta, hOpts)
        cdef HarmonicModeContainerWrapper modeWrap = HarmonicModeContainerWrapper()
        modeWrap.wrap(modescpp)

        cdef int[::1] l = modeWrap.lmodes.data
        cdef int[::1] m = modeWrap.mmodes.data
        cdef int modeNum = len(l)

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(timeSteps, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(timeSteps, dtype=np.float64)
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross)

        self.hcpp.computeWaveformHarmonics(dereference(h.hcpp), &l[0], &m[0], modeNum, dereference(inspiral.inspiralcpp), theta, phi - Phi_phi0)
        if return_list:
            return [plus, cross]

        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] waveform = np.empty(timeSteps, dtype=np.complex128)
        waveform = -1.j*cross
        waveform += plus

        return waveform