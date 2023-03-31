from bhpwaveformcy import HarmonicAmplitudesPy
import numpy as np

class HarmonicAmplitudes:
    def __init__(self, lmodes = None, mmodes = None):
        if isinstance(lmodes, np.ndarray) and isinstance(mmodes, np.ndarray):
            if lmodes.shape == mmodes.shape:
                self.amplitude_generator = HarmonicAmplitudesPy(lmodes, mmodes, dealloc_flag = False)
            else: 
                self.amplitude_generator = HarmonicAmplitudesPy(dealloc_flag = False)
        else:
            self.amplitude_generator = HarmonicAmplitudesPy(dealloc_flag = False)

    def __call__(self, l, m, a, r):
        if isinstance(a, list) or isinstance(a, np.ndarray) or isinstance(r, list) or isinstance(r, np.ndarray):
            return self.amplitude_generator.amplitude_array(l, m, a, r)*np.exp(1.j*self.amplitude_generator.phase_array(l, m, a, r))
        
        return self.amplitude_generator.amplitude(l, m, a, r)*np.exp(1.j*self.amplitude_generator.phase(l, m, a, r))