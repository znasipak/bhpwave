from libcpp.vector cimport vector
import numpy as np
cimport numpy as np
from libcpp.string cimport string
from libcpp.complex cimport complex as cpp_complex

include "harmonic_wrap.pyx"

cdef extern from "waveform.hpp":
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

    cdef cppclass WaveformHarmonicOptions:
        WaveformHarmonicOptions()
        WaveformHarmonicOptions(double rescale, int num, int pad_output)

        cpp_complex[double] rescale
        int num_threads
        int pad_output

    cdef cppclass WaveformHarmonicGenerator:
        WaveformHarmonicGenerator(HarmonicAmplitudes &Alm, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts) except +

        void computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi)
        void computeWaveformHarmonic(WaveformContainer &h, int l, int m, InspiralContainer &inspiral, double theta, double phi)
        void computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi) except +

        void computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts)
        void computeWaveformHarmonic(WaveformContainer &h, int l, int m, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts)
        void computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts) except +
        
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
        WaveformHarmonicOptions getWaveformHarmonicOptions()
        HarmonicOptions getHarmonicOptions()

cdef class WaveformContainerNumpyWrapper:
    cdef WaveformContainer *hcpp

    def __cinit__(self, np.ndarray[ndim=1, dtype=np.float64_t, mode='c'] plus, np.ndarray[ndim=1, dtype=np.float64_t, mode='c'] cross, int time_steps):
        self.hcpp = new WaveformContainer(&plus[0], &cross[0], time_steps)

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

cdef class WaveformHarmonicGeneratorPyWrapper:
    cdef WaveformHarmonicGenerator *hcpp
    cdef int modeCheck

    def __cinit__(self, HarmonicAmplitudesPy Alm, dict harmonic_kwargs = {}, dict waveform_kwargs = {}):
        cdef WaveformHarmonicOptions wOpts
        cdef HarmonicOptions hOpts

        if "check_modes" in harmonic_kwargs.keys():
            self.modeCheck = harmonic_kwargs["check_modes"]
        else:
            self.modeCheck = 0
        
        if "eps" in harmonic_kwargs.keys():
            hOpts.epsilon = harmonic_kwargs["eps"]
        if "max_samples" in harmonic_kwargs.keys():
            hOpts.max_samples = harmonic_kwargs["max_samples"]
        
        if "num_threads" in waveform_kwargs.keys():
            wOpts.num_threads = waveform_kwargs["num_threads"]
        if "pad_output" in waveform_kwargs.keys():
            wOpts.pad_output = waveform_kwargs["pad_output"]

        self.hcpp = new WaveformHarmonicGenerator(dereference(Alm.harmonicscpp), hOpts, wOpts)

    def __dealloc__(self):
        del self.hcpp
 
    cdef evaluate_harmonics(self, np.ndarray[ndim = 1, dtype = int, mode='c'] lmodes, np.ndarray[ndim = 1, dtype = int, mode='c'] mmodes, InspiralContainerWrapper inspiral, double theta, double phi):
        cdef WaveformContainerWrapper h = WaveformContainerWrapper(inspiral.timesteps)
        self.hcpp.computeWaveformHarmonics(dereference(h.hcpp), &lmodes[0], &mmodes[0], lmodes.shape[0], dereference(inspiral.inspiralcpp), theta, phi)
        return h

    def __call__(self, l, m, InspiralContainerWrapper inspiral, double theta, double phi):
        cdef WaveformContainerWrapper h = WaveformContainerWrapper(inspiral.timesteps)
        if isinstance(l, list) or isinstance(m, list):
            l = np.array(l, dtype=np.int32)
            m = np.array(m, dtype=np.int32)
        if isinstance(l, np.ndarray) and isinstance(m, np.ndarray):
            assert (l.shape == m.shape), "Shapes of {}, {} for lmodes and mmodes are incompatible".format(l.shape, m.shape)
            return self.evaluate_harmonics(l, m, inspiral, theta, phi)
        else:
            self.hcpp.computeWaveformHarmonic(dereference(h.hcpp), l, m, dereference(inspiral.inspiralcpp), theta, phi)
        return h


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

        self.hcpp = new WaveformGenerator(dereference(traj.trajcpp), dereference(Alm.harmonicscpp), hOpts, wOpts)

    def __dealloc__(self):
        del self.hcpp

    def time_step_number(self, double M, double mu, double a, double r0, double dt, double T):
        return self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)

    def select_modes(self, double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, pad_nmodes = False, **kwargs):
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

    def waveform_harmonics(self, int[::1] l, int[::1] m, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
        
        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(timeSteps, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(timeSteps, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] waveform = np.empty(timeSteps, dtype=np.complex128)
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross, timeSteps)
        
        self.hcpp.computeWaveform(dereference(h.hcpp), &l[0], &m[0], l.shape[0], M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, hOpts, wOpts)
        waveform = -1.j*cross
        waveform += plus

        return waveform
    
    def waveform(self, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, bint pad_output = False, bint return_list = False, **kwargs):
        cdef int timeSteps
        cdef WaveformHarmonicOptions wOpts = self.hcpp.getWaveformHarmonicOptions()
        cdef HarmonicOptions hOpts = self.hcpp.getHarmonicOptions()
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
        
        if "eps" in kwargs.keys():
            hOpts.epsilon = kwargs["eps"]
        if "max_samples" in kwargs.keys():
            hOpts.max_samples = kwargs["max_samples"]

        if "num_threads" in kwargs.keys():
            wOpts.num_threads = kwargs["num_threads"]
        # cdef WaveformContainerWrapper h = WaveformContainerWrapper(timeSteps)

        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] plus = np.zeros(timeSteps, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.float64_t, mode='c'] cross = np.zeros(timeSteps, dtype=np.float64)
        cdef np.ndarray[ndim = 1, dtype = np.complex128_t, mode='c'] waveform = np.empty(timeSteps, dtype=np.complex128)
        cdef WaveformContainerNumpyWrapper h = WaveformContainerNumpyWrapper(plus, cross, timeSteps)
        
        self.hcpp.computeWaveform(dereference(h.hcpp), M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, hOpts, wOpts)

        if return_list:
            return [plus, cross]
        else:
            waveform = -1.j*cross
            waveform += plus

        return waveform

    def waveform_source_frame(self, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double dt, double T, bint pad_output = False):
        cdef int timeSteps
        if pad_output:
            timeSteps = self.hcpp.computeTimeStepNumber(dt, T)
        else:
            timeSteps = self.hcpp.computeTimeStepNumber(M, mu, a, r0, dt, T)
        cdef WaveformContainerWrapper h = WaveformContainerWrapper(timeSteps)
        self.hcpp.computeWaveformSourceFrame(dereference(h.hcpp), M, mu, a, r0, theta, phi, Phi_phi0, dt, T)
        return h
