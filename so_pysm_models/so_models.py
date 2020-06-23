try:  # PySM >= 3.2.1
    import pysm3.units as u
    import pysm3 as pysm
except ImportError:
    import pysm.units as u
    import pysm

from .utils import RemoteData

def get_so_models(key, nside, map_dist=None, coord="C"):

    remote_data = RemoteData(coord)

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
            dust_T = remote_data.get("dust_T_ns{}.fits".format(nside_template))
            freq_ref_I = 545 * u.GHz
        elif key == "SO_d1":
            map_mbb_index = remote_data.get(
                "variable_spectral_index/beta_dust_ns{}_1deg.fits".format(
                    nside_template
                )
            )
            map_mbb_temperature = remote_data.get(
                "variable_spectral_index/temperature_dust_ns{}_1deg.fits".format(
                    nside_template
                )
            )
            dust_T = remote_data.get("dust_T_ns{}_353GHz.fits".format(nside_template))
            freq_ref_I = 353 * u.GHz
        model = pysm.ModifiedBlackBody(
            dust_T,
            remote_data.get("dust_Q_ns{}.fits".format(nside_template)),
            remote_data.get("dust_U_ns{}.fits".format(nside_template)),
            freq_ref_I=freq_ref_I,
            freq_ref_P=353 * u.GHz,
            map_mbb_index=map_mbb_index,
            map_mbb_temperature=map_mbb_temperature,
            nside=nside,
            has_polarization=True,
            unit_I=u.uK_RJ,
            unit_Q=u.uK_RJ,
            unit_U=u.uK_RJ,
            unit_mbb_temperature=u.K,
            map_dist=map_dist,
        )
    elif key.startswith("SO_s"):
        if key == "SO_s0":
            index = -3.1
            synch_T = remote_data.get("synch_T_ns{}.fits".format(nside_template))
            freq_ref_I = 408 * u.MHz
        elif key == "SO_s1":
            index = remote_data.get(
                "variable_spectral_index/beta_synch_ns{}_1deg.fits".format(
                    nside_template
                )
            )
            synch_T = remote_data.get("synch_T_ns{}_23GHz.fits".format(nside_template))
            freq_ref_I = 23 * u.GHz
        model = pysm.PowerLaw(
            map_I=synch_T,
            map_Q=remote_data.get("synch_Q_ns{}.fits".format(nside_template)),
            map_U=remote_data.get("synch_U_ns{}.fits".format(nside_template)),
            freq_ref_I=freq_ref_I,
            freq_ref_P=23 * u.GHz,
            map_pl_index=index,
            nside=nside,
            unit_I=u.uK_RJ,
            unit_Q=u.uK_RJ,
            unit_U=u.uK_RJ,
            map_dist=map_dist,
        )
    elif key == "SO_f0":
        model = pysm.PowerLaw(
            map_I=remote_data.get("freefree_T_ns{}.fits".format(nside_template)),
            freq_ref_I=30 * u.GHz,
            map_pl_index=-2.14,
            nside=nside,
            unit_I=u.uK_RJ,
            map_dist=map_dist,
        )
    elif key == "SO_a0":
        model = pysm.Sky(
            component_objects=[
                pysm.SpDust(
                    map_I=remote_data.get("ame1_T_ns{}.fits".format(nside_template)),
                    freq_ref_I=22.8 * u.GHz,
                    emissivity=remote_data.get("ame_emissivity.txt"),
                    freq_peak=18.95 * u.GHz,
                    freq_ref_peak=30 * u.GHz,
                    nside=nside,
                    unit_I=u.uK_RJ,
                    map_dist=map_dist,
                ),
                pysm.SpDust(
                    map_I=remote_data.get("ame2_T_ns{}.fits".format(nside_template)),
                    freq_ref_I=41.0 * u.GHz,
                    emissivity=remote_data.get("ame_emissivity.txt"),
                    freq_peak=33.35 * u.GHz,
                    freq_ref_peak=30 * u.GHz,
                    nside=nside,
                    unit_I=u.uK_RJ,
                    map_dist=map_dist,
                ),
            ]
        )
    elif key == "SO_a1":
        model = pysm.Sky(
            component_objects=[
                pysm.SpDustPol(
                    map_I=remote_data.get("ame1_T_ns{}.fits".format(nside_template)),
                    freq_ref_I=22.8 * u.GHz,
                    emissivity=remote_data.get("ame_emissivity.txt"),
                    freq_peak=remote_data.get(
                        "variable_spectral_index/ame_nu0_peak_ns{}_1deg.fits".format(
                            nside_template
                        )
                    ),
                    freq_ref_peak=30 * u.GHz,
                    pol_frac=0.01,
                    angle_Q=remote_data.get("dust_Q_ns{}.fits".format(nside_template)),
                    angle_U=remote_data.get("dust_U_ns{}.fits".format(nside_template)),
                    nside=nside,
                    unit_I=u.uK_RJ,
                    map_dist=map_dist,
                ),
                pysm.SpDustPol(
                    map_I=remote_data.get("ame2_T_ns{}.fits".format(nside_template)),
                    freq_ref_I=41.0 * u.GHz,
                    emissivity=remote_data.get("ame_emissivity.txt"),
                    freq_peak=33.35 * u.GHz,
                    freq_ref_peak=30 * u.GHz,
                    pol_frac=0.01,
                    angle_Q=remote_data.get("dust_Q_ns{}.fits".format(nside_template)),
                    angle_U=remote_data.get("dust_U_ns{}.fits".format(nside_template)),
                    nside=nside,
                    unit_I=u.uK_RJ,
                    map_dist=map_dist,
                ),
            ]
        )

    return model
