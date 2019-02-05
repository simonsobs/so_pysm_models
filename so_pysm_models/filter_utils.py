import numpy as np
import healpy as hp

def create_high_pass_filter(l1, l2, lmax):
    ell = np.arange(lmax+1)
    wl = 0.5 * (1 - np.cos(np.pi * (l2 - ell) / (l2 - l1)))
    wl[ell >= l2] = 1.0
    return wl


def create_low_pass_filter(l1, l2, lmax):
    ell = np.arange(lmax+1)
    wl = 0.5 * (1 - np.cos(np.pi * (l2 - ell) / (l2 - l1)))
    wl[ell <= l1] = 1.0
    return wl


def apply_filter(hmap, filt):
    nside = hp.get_nside(hmap)
    lmax = 3 * nside
    filt = filt[0:lmax]
    alm = hp.map2alm(hmap, lmax=lmax)
    almf = hp.almxfl(alm, filt)
    hmap_out = hp.alm2map(almf, nside)
    return hmap_out
