from bhpswshcy import Yslm as YslmCy
import numpy as np

def Yslm(s, l, m, theta, phi):
    """The spin-weighted spherical harmonics 
    :math:`{}_s Y_{lm}(\\theta, \\phi) = {}_s P_{lm}(\\cos\\theta)e^{im\\phi}`,
    where :math:`{}_s P_{lm}(z)` are the spin-weighted Legendre polynomials.

    Parameters
    ----------
    s : int
        spin-weight
    l : int
        polar mode number
    m : int
        azimuthal mode number
    theta : double or array[double]
        polar angle
    phi : double or array[double]
        azimuthal angle
    """
    return YslmCy(s, l, m, theta)*np.exp(1.j*m*phi)

def Pslm(s, l, m, z):
    """The spin-weighted Legendre function :math:`{}_s P_{lm}(z)`.

    Parameters
    ----------
    s : int
        spin-weight
    l : int
        polar mode number
    m : int
        azimuthal mode number
    z : double or array[double]
        argument :math:`-1 \\leq z \\leq 1`
    """
    theta = np.arccos(z)
    return YslmCy(s, l, m, theta)