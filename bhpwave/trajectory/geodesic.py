import numpy as np

######################################
# Define Useful Trajectory Functions #
######################################

OMEGA_MIN = 2.e-3
A_MAX = 0.9999

def kerr_circ_geo_radius(a, omega):
    """Calculates the Boyer-Lindquist radius of a circular geodesic with orbital
    frequency :code:`omega` in a Kerr spacetime paramtrized by the Kerr spin :code:`a`

    Parameters
    ----------
    a : double or array
        Kerr spin parameter
    omega : double or array
        orbital frequency

    Returns
    -------
    double or array
    """
    return (abs(omega)*(1. - a*omega)/(omega**2))**(2./3.)

def kerr_circ_geo_orbital_frequency(a, r):
    """Calculates the orbital frequency of a circular geodesic with Boyer-Lindquist radius
    :code:`r` in a Kerr spacetime paramtrized by the Kerr spin :code:`a`

    Parameters
    ----------
    a : double or array
        Kerr spin parameter
    r : double or array
        orbital radius

    Returns
    -------
    double or array
    """
    v = 1./np.sqrt(r)
    return pow(v, 3)/(1 + a*pow(v, 3))

def kerr_isco_radius(a):
    """Calculates the Boyer-Lindquist radius of the innermost stable circular orbit (ISCO)
    in a Kerr spacetime paramtrized by the Kerr spin :code:`a`

    Parameters
    ----------
    a : double or array
        Kerr spin parameter

    Returns
    -------
    double or array
    """
    sgnX = np.sign(a)
    z1 = 1 + pow(1 - a*a, 1./3.)*(pow(1 - a, 1./3.) + pow(1 + a, 1./3.))
    z2 = np.sqrt(3*a*a + z1*z1)

    return 3 + z2 - sgnX*np.sqrt((3. - z1)*(3. + z1 + 2.*z2))

def kerr_isco_frequency(a):
    """Calculates the orbital frequency of the innermost stable circular orbit (ISCO)
    in a Kerr spacetime paramtrized by the Kerr spin :code:`a`

    Parameters
    ----------
    a : double or array
        Kerr spin parameter

    Returns
    -------
    double or array
    """
    rISCO = kerr_isco_radius(a)
    return kerr_circ_geo_orbital_frequency(a, rISCO)

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