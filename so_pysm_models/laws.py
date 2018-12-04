import numpy as np


def power_law(nu, nu_ref_s, beta_s):
    """ Function to compute power law SED.
    Parameters
    ----------
    nu: float, or array_like(float)
        Frequency in GHz.
    beta_s: float
        Power law index in RJ units.
    Returns
    -------
    array_like(float)
        SED relative to reference frequency.
    """
    x = nu / nu_ref_s
    sed = x ** beta_s
    return sed


def curved_power_law(nu, nu_ref_s, beta_s, beta_c):
    """ Function to compute curved power law SED.
    Parameters
    ----------
    nu: float, or array_like(float)
        Frequency in GHz.
    beta_s: float
        Power law index in RJ units.
    beta_c: float
        Power law index curvature.
    Returns
    -------
    array_like(float)
        SED relative to reference frequency.
    """
    x = nu / nu_ref_s
    sed = x ** (beta_s + beta_c * np.log(nu / nu_ref_s))
    return sed


def black_body_cmb(nu):
    """ Function to compute CMB SED.
    Parameters
    ----------
    nu: float, or array_like(float)
        Frequency in GHz.
    """
    x = 0.0176086761 * nu
    ex = np.exp(x)
    sed = ex * (x / (ex - 1)) ** 2
    return sed
