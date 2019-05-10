import healpy as hp
import numpy as np
import pysm
import os
from astropy.utils import data
from pysm import read_map
import pysm.units as u
def loadtxt(*args):
    pass

DATAURL = "http://portal.nersc.gov/project/cmb/so_pysm_models_data/"
PREDEFINED_DATA_FOLDERS = ["/global/project/projectdirs/cmb/www/so_pysm_models_data"]


def get_data_from_url(filename):
    """Retrieves input templates from remote server,
    in case data is available in one of the PREDEFINED_DATA_FOLDERS defined above,
    e.g. at NERSC, those are directly returned."""
    for folder in PREDEFINED_DATA_FOLDERS:
        full_path = os.path.join(folder, filename)
        if os.path.exists(full_path):
            return full_path
    with data.conf.set_temp("dataurl", DATAURL), data.conf.set_temp(
        "remote_timeout", 30
    ):
        map_out = data.get_pkg_data_filename(filename, show_progress=True)
    return map_out


def get_so_models(key, nside, pixel_indices=None, mpi_comm=None):
    small_scale = key.endswith("s")
    if small_scale:
        nside_template = 4096
        key = key[:-1]
    else:
        nside_template = 512

    if key == "SO_d0":
        model = pysm.ModifiedBlackBody(
            get_data_from_url("dust_T_ns{}.fits".format(nside_template)),
            get_data_from_url("dust_Q_ns{}.fits".format(nside_template)),
            get_data_from_url("dust_U_ns{}.fits".format(nside_template)),
            freq_ref_I=545 * u.GHz,
            freq_ref_P=353 * u.GHz,
            map_mbb_index=1.53,
            map_mbb_temperature=19.6 * u.K,
            nside=nside,
            has_polarization=True,
            unit_I=u.uK_RJ,
            unit_Q=u.uK_RJ,
            unit_U=u.uK_RJ,
            unit_mbb_temperature=u.K,
            pixel_indices=pixel_indices,
            mpi_comm=mpi_comm,
        )
    return model


def SO_d1(nside, pixel_indices=None, mpi_comm=None, nside_template=512):
    T_map = get_data_from_url("dust_T_ns{}.fits".format(nside_template))
    Q_map = get_data_from_url("dust_Q_ns{}.fits".format(nside_template))
    U_map = get_data_from_url("dust_U_ns{}.fits".format(nside_template))
    beta_map = get_data_from_url("variable_spectral_index/beta_dust_ns{}_1deg.fits".format(nside_template))
    temperature_map = get_data_from_url("variable_spectral_index/temperature_dust_ns{}_1deg.fits".format(nside_template))
    A_I = read_map(
        T_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
    )
    return [
        {
            "model": "modified_black_body",
            "nu_0_I": 545.,
            "nu_0_P": 353.,
            "A_I": A_I,
            "A_Q": read_map(
                Q_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "A_U": read_map(
                U_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "spectral_index": read_map(
                beta_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "temp": read_map(
                temperature_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "add_decorrelation": False,
        }
    ]


def SO_s0(nside, pixel_indices=None, mpi_comm=None, nside_template=512):
    T_map = get_data_from_url("synch_T_ns{}.fits".format(nside_template))
    Q_map = get_data_from_url("synch_Q_ns{}.fits".format(nside_template))
    U_map = get_data_from_url("synch_U_ns{}.fits".format(nside_template))
    A_I = read_map(
        T_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
    )
    return [
        {
            "model": "power_law",
            "nu_0_I": 0.408,
            "nu_0_P": 23.,
            "A_I": A_I,
            "A_Q": read_map(
                Q_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "A_U": read_map(
                U_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "spectral_index": np.ones(len(A_I)) * -3.1,
        }
    ]

def SO_s1(nside, pixel_indices=None, mpi_comm=None, nside_template=512):
    T_map = get_data_from_url("synch_T_ns{}.fits".format(nside_template))
    Q_map = get_data_from_url("synch_Q_ns{}.fits".format(nside_template))
    U_map = get_data_from_url("synch_U_ns{}.fits".format(nside_template))
    beta_map = get_data_from_url("variable_spectral_index/beta_synch_ns{}_1deg.fits".format(nside_template))
    A_I = read_map(
        T_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
    )
    return [
        {
            "model": "power_law",
            "nu_0_I": 0.408,
            "nu_0_P": 23.,
            "A_I": A_I,
            "A_Q": read_map(
                Q_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "A_U": read_map(
                U_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "spectral_index": read_map(
                beta_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
        }
    ]

def SO_f0(nside, pixel_indices=None, mpi_comm=None, nside_template=512):
    T_map = get_data_from_url("freefree_T_ns{}.fits".format(nside_template))
    return [
        {
            "model": "power_law",
            "nu_0_I": 30.,
            "A_I": read_map(
                T_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "spectral_index": -2.14,
        }
    ]


def SO_a0(nside, pixel_indices=None, mpi_comm=None, nside_template=512):
    T_map1 = get_data_from_url("ame1_T_ns{}.fits".format(nside_template))
    T_map2 = get_data_from_url("ame2_T_ns{}.fits".format(nside_template))
    return [
        {
            "model": "spdust",
            "nu_0_I": 22.8,
            "nu_0_P": 22.8,
            "A_I": read_map(
                T_map1, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "nu_peak_0": 30.,
            "emissivity": loadtxt(
                get_data_from_url("ame_emissivity.txt"), mpi_comm=mpi_comm, unpack=True
            ),
            "nu_peak": 18.95,
        },
        {
            "model": "spdust",
            "nu_0_I": 41.0,
            "nu_0_P": 41.0,
            "A_I": read_map(
                T_map2, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "nu_peak_0": 30.,
            "emissivity": loadtxt(
                get_data_from_url("ame_emissivity.txt"), mpi_comm=mpi_comm, unpack=True
            ),
            "nu_peak": 33.35,
        },
    ]

def SO_a1(nside, pixel_indices=None, mpi_comm=None, nside_template=512):
    T_map1 = get_data_from_url("ame1_T_ns{}.fits".format(nside_template))
    T_map2 = get_data_from_url("ame2_T_ns{}.fits".format(nside_template))
    Q_dust = get_data_from_url("dust_Q_ns{}.fits".format(nside_template))
    U_dust = get_data_from_url("dust_U_ns{}.fits".format(nside_template))
    nu_peak_map = get_data_from_url("variable_spectral_index/ame_nu0_peak_ns{}_1deg.fits".format(nside_template))
    return [
        {
            "model": "spdust_pol",
            "nu_0_I": 22.8,
            "nu_0_P": 22.8,
            "A_I": read_map(
                T_map1, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "nu_peak_0": 30.,
            "emissivity": loadtxt(
                get_data_from_url("ame_emissivity.txt"), mpi_comm=mpi_comm, unpack=True
            ),
            'nu_peak': read_map(
                nu_peak_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            'pol_frac': 0.01,
            'angle_q': read_map(
                Q_dust, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            'angle_u': read_map(
                U_dust, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            )
        },
        {
            "model": "spdust_pol",
            "nu_0_I": 41.0,
            "nu_0_P": 41.0,
            "A_I": read_map(
                T_map2, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            "nu_peak_0": 30.,
            "emissivity": loadtxt(
                get_data_from_url("ame_emissivity.txt"), mpi_comm=mpi_comm, unpack=True
            ),
            "nu_peak": 33.35,
            'pol_frac': 0.01,
            'angle_q': read_map(
                Q_dust, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            ),
            'angle_u': read_map(
                U_dust, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm
            )
        },
    ]
