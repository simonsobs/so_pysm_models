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


def modified_black_body(nu, nu_ref_d, beta_d, T_d):
    """ Function to compute modified blackbody dust SED.
    Parameters
    ----------
    nu: float or array_like(float)
        Frequency at which to calculate SED.
    nu_ref_d: float
        Reference frequency in GHz.
    beta_d: float
        Power law index of dust opacity.
    T_d: float
        Temperature of the dust.
    Returns
    -------
    array_like(float)
        SED of dust modified black body relative to reference frequency.
    """
    x_to = 0.0479924466 * nu / T_d
    x_from = 0.0479924466 * nu_ref_d / T_d
    sed = (nu / nu_ref_d) ** (1 + beta_d) * (np.exp(x_from) - 1) / (np.exp(x_to) - 1)
    return sed
