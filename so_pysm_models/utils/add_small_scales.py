import healpy as hp
import numpy as np

def create_high_pass_filter(l1, l2, lmax):
    ell = np.arange(l1, l2)
    wl = np.zeros(lmax+1)
    wl[l2:] = 1.0
    wl[ell] = 0.5 * (1 - np.cos(np.pi * (ell - l1) / (l2 - l1)))
    return wl


def create_low_pass_filter(l1, l2, lmax):
    ell = np.arange(l1, l2)
    wl = np.zeros(lmax+1)
    wl[0:l1] = 1.0
    wl[ell] = 0.5 * (1 - np.cos(np.pi * (l2 - ell) / (l2 - l1)))
    return wl

def apply_filter(hmap, filt):
    import healpy as hp
    nside = hp.get_nside(hmap)
    print('apply_filter', nside)
    lmax = 3*nside
    filt = filt[0:lmax]
    alm = hp.map2alm(hmap, lmax=lmax)
    almf = np.array([hp.almxfl(alm[i], filt) for i in range(len(alm))])
    hmap_out = hp.alm2map(almf, nside)
    return hmap_out

def add_gaussian_small_scales(map_in, nside_out, pol=False):
    cl_in = hp.anafast(map_in)
    map_in_ns_out = hp.ud_grade(map_in, nside_out)
    map_in_ns_out_smt = hp.smoothing(map_in_ns_out, fwhm=np.radians(180./1000.))
    print('map smoothed')
    ell_to_fit = np.arange(100, 500)
    hpf = create_high_pass_filter(300, 1000, 4*nside_out)
    ell_hell = np.arange(len(hpf))
    if pol==False:
    	cl_T_to_fit  = cl_in[ell_to_fit]
    else:
    	cl_T_to_fit = cl_in[0][ell_to_fit]
    fit_cl_T = np.polyfit(np.log(ell_to_fit), np.log(cl_T_to_fit), 1)
    cl_hell_T = np.exp(fit_cl_T[0]*np.log(ell_hell)+fit_cl_T[1])*hpf
    cl_hell_T[0] = 0
    cl_hell_zero = np.zeros(len(cl_hell_T))
    if pol:
        cl_E_to_fit = cl_in[1][ell_to_fit]
        cl_B_to_fit = cl_in[2][ell_to_fit]
        fit_cl_E = np.polyfit(np.log(ell_to_fit), np.log(cl_E_to_fit), 1)
        fit_cl_B = np.polyfit(np.log(ell_to_fit), np.log(cl_B_to_fit), 1)
        cl_hell_E = np.exp(fit_cl_E[0]*np.log(ell_hell)+fit_cl_E[1])*hpf
        cl_hell_B = np.exp(fit_cl_B[0]*np.log(ell_hell)+fit_cl_B[1])*hpf
        cl_hell_E[0] = 0
        cl_hell_B[0] = 0
        cl_hell = np.array([cl_hell_T, cl_hell_E, cl_hell_B, cl_hell_zero, cl_hell_zero, cl_hell_zero])
        map_ss = hp.synfast(cl_hell, nside_out, lmax=nside_out*3, pol=True, new=True)
    else:
        map_ss = hp.synfast(cl_hell_T, nside_out, lmax=nside_out*3)
    print('map small scales computed')
    map_ss_mod = map_ss*map_in_ns_out_smt
    if pol==False:
    	coeff_T = np.std(map_ss_mod)/np.std(map_ss)
    	map_out_T = map_ss_mod/coeff_T+map_in_ns_out_smt
    else:
    	coeff_T = np.std(map_ss_mod[0])/np.std(map_ss[0])
    	map_out_T = map_ss_mod[0]/coeff_T+map_in_ns_out_smt[0]
    if np.any(map_out_T<0):
        negative_pix = np.where(map_out_T<0)[0]
        print('negative pixels ', len(negative_pix))
        if pol==False:
            map_out_T[negative_pix] = map_in_ns_out[negative_pix]
        else:
            map_out_T[negative_pix] = map_in_ns_out[0][negative_pix]
    if pol:
        coeff_Q = np.std(map_ss_mod[1])/np.std(map_ss[1])
        coeff_U = np.std(map_ss_mod[2])/np.std(map_ss[2])
        coeff_P = (coeff_Q+coeff_U)/2.
        map_out_Q = map_ss_mod[1]/coeff_P+map_in_ns_out_smt[1]
        map_out_U = map_ss_mod[2]/coeff_P+map_in_ns_out_smt[2]
        map_out = np.array([map_out_T, map_out_Q, map_out_U])
    else:
        map_out = map_out_T
    return map_out
