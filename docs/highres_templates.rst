High resolution templates
*************************

Starting from version 2.0, all input templates have been rotated to Equatorial, therefore by default the output is in Equatorial coordinates.

:py:mod:`so_pysm_models` also provides access to templates with higher resolution and with updated
data compared to the models included in PySM.

They can be accessed with the function :py:func:`.get_so_models` which works similarly to the `models`
function available in PySM, for example::

    from so_pysm_models import get_so_models
    from pysm.nominal import models
    from pysm import Sky
    sky = Sky(component_objects=[get_so_models("SO_d0s", nside=4096)])

Consider that the :py:mod:`so_pysm_models` models are by default in **Equatorial coordinates**, therefore they should not
be mixed in the same run with standard PySM components which instead are in **Galactic coordinates**. If you need PySM to simultaneously handle inputs in different reference frames, please `open an issue in the PySM repository <https://github.com/healpy/pysm/issues/>`_.

:py:mod:`so_pysm_models` retrieves the templates when needed from NERSC via web accessing:
http://portal.nersc.gov/project/cmb/so_pysm_models_data/
Downloaded files are stored in the `astropy` cache, generally `~/.astropy/cache` and are accessible using :py:mod:`astropy.utils.data`, e.g. :py:func:`astropy.utils.data.get_cached_urls` gives the list of downloaded files. If running at NERSC, the module automatically uses the files accessible locally from the `/project` filesystem.

Low-resolution templates are standard PySM ones at :math:`N_{side}` 512, often with updated parameters based on Planck results.
High-resolution templates are computed from the low-resolution ones, by extrapolating
power spectra considering a simple power law model, and by generating small scales as Gaussian realization of these spectra.
High-resolution templates therefore have Gaussian small scales (for :math:`\ell > ~ 1000`) modulated with large scale signal
for both temperature and polarization.

You can access the high resolution parameters at :math:`N_{side}` 4096 appending `s` (for small scale) at the end of each model name, for example::

    from so_pysm_models import get_so_models
    from pysm import Sky
    sky_highres = Sky(component_objects=[get_so_models("SO_d0s", nside=4096)])

Whatever the :math:`N_{side}` of the input model and the requested :math:`N_{side}` in :py:func:`.get_so_models`, PySM will automatically use :py:func:`healpy.ud_grade` to adjust the map resolution.


Details about individual models
===============================

Append "s" after a model name to access the :math:`N_{side}` 4096 template, i.e. `SO_f0s`.

**Dust**

* **SO_d0**: Thermal dust is modeled as a single-component modified black body, with same templates as in PySM model `d1`.  There is no spatial variation of temperature and emissivity in the sky: :math:`T=19.6` K and :math:`\beta_d=1.53` (values taken from Planck Collaboration IX 2018).

* **SO_d1**: Thermal dust is modeled as a single-component modified black body, with same templates as in PySM model `d1`.  Both spectral index and dust temperature are spatially varying up to the degree scale.

**Synchrotron**

* **SO_s0**: Templates from PySM model `s1`. Power law spectral energy distribution, with fixed spectral index :math:`\beta_s=-3.1` (from Planck Collaboration IX 2018).

* **SO_s1**: Templates from PySM model `s1`. Power law spectral energy distribution, with spatially varying spectral index up to the degree scale.

**Free Free**

* **SO_f0**: same model as PySM `f1`, no spatial variation of spectral index equal to -2.4.

**AME**

* **SO_a0**: sum of two spinning dust populations (as in PySM model `a1`) with spatially constant peak frequency. No polarization.

* **SO_a1**: sum of two spinning dust populations (as in PySM model `a1`). First one with spatially constant peak frequency, the other with spatially variable peak frequency up to the degree scale. Polarized maps simulated with thermal dust angles and nominal AME intensity, scaled globally by 1% polarization fraction.
