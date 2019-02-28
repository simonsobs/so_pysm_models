import healpy as hp
import numpy as np
import pysm
import os
from astropy.utils import data
from pysm.common import read_map, loadtxt

data_dir = os.path.join(os.path.dirname(__file__), 'so_template')
template = lambda x: os.path.join(data_dir, x)
DATAURL = "http://portal.nersc.gov/project/cmb/so_pysm_models_data/"

def get_data_from_url(filename):
    with data.conf.set_temp("dataurl", DATAURL), data.conf.set_temp(
        "remote_timeout", 30):
        map_out = data.get_pkg_data_filename(filename)
    return map_out

def models(key, nside, pixel_indices=None, mpi_comm=None, small_scale=False):
     model = eval(key)(nside, pixel_indices=pixel_indices, mpi_comm=mpi_comm, small_scale=small_scale)
     for m in model:
         m['pixel_indices'] = pixel_indices # include pixel indices in the model dictionary
         m['nside'] = nside
     return model

def SO_d0(nside, pixel_indices=None, mpi_comm=None, small_scale=False):
    if small_scale:
        T_map = get_data_from_url('dust_T_ns4096.fits')
        Q_map = get_data_from_url('dust_Q_ns4096.fits')
        U_map = get_data_from_url('dust_U_ns4096.fits')
    else:
        T_map = template('dust_T_ns512.fits')
        Q_map = template('dust_Q_ns512.fits')
        U_map = template('dust_U_ns512.fits')
    A_I = read_map(T_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm)
    return [{
        'model': 'modified_black_body',
        'nu_0_I': 545.,
        'nu_0_P': 353.,
        'A_I': A_I,
        'A_Q': read_map(Q_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'A_U': read_map(U_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'spectral_index': np.ones(len(A_I)) * 1.53,
        'temp': np.ones(len(A_I)) * 19.6,
        'add_decorrelation': False,
    }]

def SO_s0(nside, pixel_indices=None, mpi_comm=None, small_scale=False):
    if small_scale:
        T_map = get_data_from_url('synch_T_ns4096.fits')
        Q_map = get_data_from_url('synch_Q_ns4096.fits')
        U_map = get_data_from_url('synch_U_ns4096.fits')
    else:
        T_map = template('synch_T_ns512.fits')
        Q_map = template('synch_Q_ns512.fits')
        U_map = template('synch_U_ns512.fits')
    A_I = read_map(T_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm)
    return [{
        'model': 'power_law',
        'nu_0_I': 0.408,
        'nu_0_P': 23.,
        'A_I': A_I,
        'A_Q': read_map(Q_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'A_U': read_map(U_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'spectral_index': np.ones(len(A_I)) * -3.1,
    }]

def SO_f0(nside, pixel_indices=None, mpi_comm=None, small_scale=False):
    if small_scale:
        T_map = get_data_from_url('freefree_T_ns4096.fits')
    else:
        T_map = template('freefree_T_ns512.fits')
    return [{
        'model': 'power_law',
        'nu_0_I': 30.,
        'A_I': read_map(T_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'spectral_index': -2.14,
    }]

def SO_a0(nside, pixel_indices=None, mpi_comm=None, small_scale=False):
    if small_scale:
        T_map1 = get_data_from_url('ame_T_ns4096.fits')
        T_map2 = get_data_from_url('ame2_T_ns4096.fits')
    else:
        T_map1 = template('ame_T_ns512.fits')
        T_map2 = template('ame2_T_ns512.fits')
    return [{
        'model': 'spdust',
        'nu_0_I': 22.8,
        'nu_0_P': 22.8,
        'A_I': read_map(T_map1, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'nu_peak_0': 30.,
        'emissivity': loadtxt(SO_template('ame_emissivity.txt'), mpi_comm=mpi_comm, unpack=True),
        'nu_peak': 18.95,
    }, {
        'model': 'spdust',
        'nu_0_I': 41.0,
        'nu_0_P': 41.0,
        'A_I': read_map(T_map2, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'nu_peak_0': 30.,
        'emissivity': loadtxt(SO_template('ame_emissivity.txt'), mpi_comm=mpi_comm, unpack=True),
        'nu_peak': 33.35,
    }]
