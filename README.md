# Interpolating component for PySM

Adds a custom emission to the sky simulated by [PySM](https://github.com/bthorne93/PySM_public)
defined as a set of template maps at pre-defined frequencies to be interpolated at the frequencies
requested through PySM.

## How to install

Checkout this repository in the folder of your code:

    git clone https://github.com/zonca/pysm_interpolating_component interpolating_component
    
From your code:

    from interpolating_component import InterpolatingComponent
    
## Inputs

A folder of maps named with their frequency in GHz with the flux in any unit supported
by PySM (e.g. `Jysr`, `MJsr`, `uK_RJ`, `K_CMB`).
They don't need to be equally spaced

For example:

    ls `cib_precomputed_maps/`
    0010.0.fits 0015.0.fits 0018.0.fits 
    
## Usage

Instantiate `InterpolatingComponent` and point it to the folder, define the unit and the target nside (same used by PySM).
It supports all `interpolation_kind` of
[`scipy.interpolate.interp1d`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html),
e.g. "nearest", "linear", "quadratic", "cubic".

    cib = InterpolatingComponent(path="cib_precomputed_maps", input_units="MJysr", target_nside=nside, interpolation_kind="linear",
                             has_polarization=False, verbose=True)
                             
Finally add the component to a PySM sky after the other components:

```
import pysm
nside = 64
sky_config = {
    'synchrotron' : models("s1", nside),
    'dust' : models("d1", nside),
    'freefree' : models("f1", nside),
    'cmb' : models("c1", nside),
    'ame' : models("a1", nside),
}
sky = pysm.Sky(sky_config)
sky.add_component("cib", cib)
```

Then you can access it as all other components:

    sky.cib(21.) # compute emission at 21 GHz in uK_RJ
    
it is also included in the total emission when you observe the PySM sky with an instrument, optionally with beam and bandpass    
