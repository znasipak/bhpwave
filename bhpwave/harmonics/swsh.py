from bhpswshcy import Yslm as YslmCy
import numpy as np

def Yslm(s, l, m, theta, phi):
    """
    The spin-weighted spherical harmonics

    :param s: spin-weight
    :type s: int
    :param l: polar mode number
    :type l: int
    :param m: azimuthal mode number
    :type m: int
    :param theta: polar angle 
    :type theta: double or array[double]
    :param phi: azimuthal angle
    :type phi: double or array[double]
    """
    return YslmCy(s, l, m, theta)*np.exp(1.j*m*phi)

def Pslm(s, l, m, z):
    """
    The spin-weighted Legendre function

    :param s: spin-weight
    :type s: int
    :param l: polar mode number
    :type l: int
    :param m: azimuthal mode number
    :type m: int
    :param z: argument :math:`-1 \\leq z \\leq 1`
    :type z: double or array[double]
    """
    theta = np.arccos(z)
    return YslmCy(s, l, m, theta)