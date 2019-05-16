import healpy as hp
import numpy as np
import pysm
import os
from astropy.utils import data
from pysm import read_map
import pysm.units as u
from .utils import get_data_from_url

def get_so_models(key, nside, pixel_indices=None, mpi_comm=None):
    small_scale = key.endswith("s")
    if small_scale:
        nside_template = 4096
        key = key[:-1]
    else:
        nside_template = 512

    if key.startswith("SO_d"):
        if key == "SO_d0":
            map_mbb_index = 1.53
            map_mbb_temperature = 19.6 * u.K
        elif key == "SO_d1":
            map_mbb_index = get_data_from_url(
                "variable_spectral_index/beta_dust_ns{}_1deg.fits".format(
                    nside_template
                )
            )
            map_mbb_temperature = get_data_from_url(
                "variable_spectral_index/temperature_dust_ns{}_1deg.fits".format(
                    nside_template
                )
            )
        model = pysm.ModifiedBlackBody(
            get_data_from_url("dust_T_ns{}.fits".format(nside_template)),
            get_data_from_url("dust_Q_ns{}.fits".format(nside_template)),
            get_data_from_url("dust_U_ns{}.fits".format(nside_template)),
            freq_ref_I=545 * u.GHz,
            freq_ref_P=353 * u.GHz,
            map_mbb_index=map_mbb_index,
            map_mbb_temperature=map_mbb_temperature,
            nside=nside,
            has_polarization=True,
            unit_I=u.uK_RJ,
            unit_Q=u.uK_RJ,
            unit_U=u.uK_RJ,
            unit_mbb_temperature=u.K,
            pixel_indices=pixel_indices,
            mpi_comm=mpi_comm,
        )
    elif key.startswith("SO_s"):
        if key == "SO_s0":
            index = -3.1
        elif key == "SO_s1":
            index = get_data_from_url(
                "variable_spectral_index/beta_synch_ns{}_1deg.fits".format(
                    nside_template
                )
            )
        model = pysm.PowerLaw(
            map_I=get_data_from_url("synch_T_ns{}.fits".format(nside_template)),
            map_Q=get_data_from_url("synch_Q_ns{}.fits".format(nside_template)),
            map_U=get_data_from_url("synch_U_ns{}.fits".format(nside_template)),
            freq_ref_I=408 * u.MHz,
            freq_ref_P=23 * u.GHz,
            map_pl_index=index,
            nside=nside,
            unit_I=u.uK_RJ,
            unit_Q=u.uK_RJ,
            unit_U=u.uK_RJ,
            pixel_indices=pixel_indices,
            mpi_comm=mpi_comm,
        )
    elif key == "SO_f0":
        model = pysm.PowerLaw(
            map_I=get_data_from_url("freefree_T_ns{}.fits".format(nside_template)),
            freq_ref_I=30 * u.GHz,
            map_pl_index=-2.14,
            nside=nside,
            unit_I=u.uK_RJ,
            pixel_indices=pixel_indices,
            mpi_comm=mpi_comm,
        )
    elif key == "SO_a0":
        model = pysm.Sky(
            component_objects=[
                pysm.SpDust(
                    map_I=get_data_from_url("ame1_T_ns{}.fits".format(nside_template)),
                    freq_ref_I=22.8 * u.GHz,
                    emissivity=get_data_from_url("ame_emissivity.txt"),
                    freq_peak=18.95 * u.GHz,
                    freq_ref_peak=30 * u.GHz,
                    nside=nside,
                    unit_I=u.uK_RJ,
                    pixel_indices=pixel_indices,
                    mpi_comm=mpi_comm,
                ),
                pysm.SpDust(
                    map_I=get_data_from_url("ame2_T_ns{}.fits".format(nside_template)),
                    freq_ref_I=41.0 * u.GHz,
                    emissivity=get_data_from_url("ame_emissivity.txt"),
                    freq_peak=33.35 * u.GHz,
                    freq_ref_peak=30 * u.GHz,
                    nside=nside,
                    unit_I=u.uK_RJ,
                    pixel_indices=pixel_indices,
                    mpi_comm=mpi_comm,
                ),
            ]
        )
    elif key == "SO_a1":
        model = pysm.Sky(
            component_objects=[
                pysm.SpDustPol(
                    map_I=get_data_from_url("ame1_T_ns{}.fits".format(nside_template)),
                    freq_ref_I=22.8 * u.GHz,
                    emissivity=get_data_from_url("ame_emissivity.txt"),
                    freq_peak=get_data_from_url(
                        "variable_spectral_index/ame_nu0_peak_ns{}_1deg.fits".format(
                            nside_template
                        )
                    ),
                    freq_ref_peak=30 * u.GHz,
                    pol_frac=0.01,
                    angle_Q=get_data_from_url(
                        "dust_Q_ns{}.fits".format(nside_template)
                    ),
                    angle_U=get_data_from_url(
                        "dust_U_ns{}.fits".format(nside_template)
                    ),
                    nside=nside,
                    unit_I=u.uK_RJ,
                    pixel_indices=pixel_indices,
                    mpi_comm=mpi_comm,
                ),
                pysm.SpDustPol(
                    map_I=get_data_from_url("ame2_T_ns{}.fits".format(nside_template)),
                    freq_ref_I=41.0 * u.GHz,
                    emissivity=get_data_from_url("ame_emissivity.txt"),
                    freq_peak=33.35 * u.GHz,
                    freq_ref_peak=30 * u.GHz,
                    pol_frac=0.01,
                    angle_Q=get_data_from_url(
                        "dust_Q_ns{}.fits".format(nside_template)
                    ),
                    angle_U=get_data_from_url(
                        "dust_U_ns{}.fits".format(nside_template)
                    ),
                    nside=nside,
                    unit_I=u.uK_RJ,
                    pixel_indices=pixel_indices,
                    mpi_comm=mpi_comm,
                ),
            ]
        )

    return model
