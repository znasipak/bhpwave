from bhpwaveformcy import HarmonicAmplitudesPy
import numpy as np
import os

path_to_file = os.path.dirname(os.path.abspath(__file__))
amp_base = path_to_file + "/../data/circ_data"

class HarmonicAmplitudes:
    """
    Class for generating the complex amplitudes :math:`H_{lm}` for a waveform produced by a particle 
    on a circular geodesic with radius :math:`r` about a Kerr black hole with spin :math:`a`.
    Amplitudes are generated for all geodesics in Kerr spacetimes :math:`|a| \\leq 0.995` 
    with orbital radii :math:`r_\\mathrm{ISCO} \\leq r \\lesssim 50 M`. Harmonic mode amplitudes can 
    also be evaluated for :math:`(l, m)` mode values up to :math:`l_{max} = 15`.

    :param lmax: The maximum l-mode to pre-compile in the amplitude generator. The value of lmax is restricted to lmax :math:`\\leq 15`
    :type lmax: int, optional
    :param select_modes: An option to specify which :math:`(l,m)` modes are evaluated by the amplitude generator.
        The default setting loads all modes in the range :math:`2 \\leq l \\leq 15`, :math:`1 \\leq |m| \\leq l`
    :type select_modes: list or array[tuple(double, double)], optional
    """
    def __init__(self, lmax = 15, select_modes = None):
        if lmax > 15:
            lmax = 15
        if select_modes is not None:
            if isinstance(select_modes, np.ndarray) or isinstance(select_modes, list):
                lmodes = []
                mmodes = []
                for mode in select_modes:
                    lmodes.append(mode[0])
                    mmodes.append(mode[1])
                lmodes = np.array(lmodes)
                mmodes = np.array(mmodes)
                modesLessThanLMAX = lmodes < lmax
                l = np.ascontiguousarray(lmodes[modesLessThanLMAX])
                m = np.ascontiguousarray(mmodes[modesLessThanLMAX])
                self.amplitude_generator = HarmonicAmplitudesPy(l, m, filebase=amp_base, dealloc_flag = False)
            else: 
                self.amplitude_generator = HarmonicAmplitudesPy(filebase=amp_base, dealloc_flag = False)
        else:
            self.amplitude_generator = HarmonicAmplitudesPy(filebase=amp_base, dealloc_flag = False)

    def __call__(self, l, m, a, r):
        """
        Generates the complex amplitude :math:`H_{lm}(a, r)`

        :param l: The spherical harmonic polar mode number
        :type l: int
        :param m: The spherical harmonic azimuthal mode number
        :type m: int
        :param a: Kerr spin parameter
        :type a: double
        :param r: Boyer-Lndquist radius of a circular geodesic in units of mass :math:`M`
        :type r: double

        :rtype: double
        """
        if isinstance(a, list) or isinstance(a, np.ndarray) or isinstance(r, list) or isinstance(r, np.ndarray):
            return self.amplitude_generator.amplitude_array(l, m, a, r)*np.exp(1.j*self.amplitude_generator.phase_array(l, m, a, r))
        
        return self.amplitude_generator.amplitude(l, m, a, r)*np.exp(1.j*self.amplitude_generator.phase(l, m, a, r))
    
    @property
    def base_class(self):
        """
        Returns the base Cython class
        """
        return self.amplitude_generator