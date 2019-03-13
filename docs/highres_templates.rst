High resolution templates
*************************

:py:mod:`so_pysm_models` also provides access to templates with higher resolution and with updated
data compared to the models included in PySM.

They can be accessed with the function :py:func:`get_so_models` which works similarly to the `models`
function available in PySM.

However it retrieves the templates from NERSC via web accessing:
http://portal.nersc.gov/project/cmb/so_pysm_models_data/

Each template has a high-resolution (:math:`N_{side}` 4096) and a low-resolution (:math:`N_{side}` 512) which are selected based
on the `small_scale` keyword passed to :py:func:`get_so_models`, i.e. you can run at :math:`N_{side}` 4096 but if you
are not interested in small scale structure you can save memory and time setting `small_scale` to `False` and
retrieve the :math:`N_{side}` 512 template.

Downloaded files are stored in the `astropy` cache, generally `~/.astropy/cache` and are accessible using :py:mod:`astropy.utils.data`, e.g. :py:func:`astropy.utils.data.get_cached_urls` gives the list of downloaded files. If running at NERSC, the module automatically uses the files accessible locally from the `/project` filesystem.


Details about individual models
===============================
