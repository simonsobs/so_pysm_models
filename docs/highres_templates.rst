High resolution templates
*************************

:py:mod:`so_pysm_models` also provides access to templates with higher resolution and with updated
data compared to the models included in PySM.

They can be accessed with the function :py:func:`get_so_models` which works similarly to the `models`
function available in PySM.

However it retrieves the templates from NERSC via web accessing:
http://portal.nersc.gov/project/cmb/so_pysm_models_data/

Downloaded files are stored in the `astropy` cache, generally `~/.astropy/cache` and are accessible using :py:mod:`astropy.utils.data`, e.g. :py:func:`astropy.utils.data.get_cached_urls` gives the list of downloaded files. If running at NERSC, the module automatically uses the files accessible locally from the `/project` filesystem.

Each template has a high-resolution (:math:`N_{side}` 4096) and a low-resolution (:math:`N_{side}` 512) which are selected based
on the `small_scale` keyword passed to :py:func:`get_so_models`, i.e. you can run at :math:`N_{side}` 4096 but if you
are not interested in small scale structure you can save memory and time setting `small_scale` to `False` and
retrieve the :math:`N_{side}` 512 template.

Low-resolution templates are standard PySM ones at :math:`N_{side}` 512.
High-resolution templates are computed from the low-resolution ones, by extrapolating
power spectra considering a simple power law model, and by generating small scales as Gaussian realization of these spectra.
High-resolution templates therefore have Gaussian small scales (for :math:`\ell > ~ 1000`) modulated with large scale signal
for both temperature and polarization.


Details about individual models
===============================

## Dust

* **SO_d0**: Thermal dust is modeled as a single-component modified black body (mbb).  There is no spatial variation of temperature and emissivity in the sky: :math:`T=19.6` K and :math:`\beta_d=1.53`.

## Synchrotron

* **SO_s0**: power law spectral energy distribution with fixed spectral index :math:`\beta_s=-3.1`.

## Free Free

* **SO_f0**: same model as PySM `f1`, no spatial variation of spectral index equal to -2.4.

## AME

* **SO_a0**: sum of two spinning dust populations with spatially constant peak frequency. No polarization.
