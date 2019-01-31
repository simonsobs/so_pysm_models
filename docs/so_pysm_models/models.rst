Summary of Models
**********************

GaussianSynchrotron
===========

This class implements Gaussian simulations for Galactic synchrotron emission.
The inputs are a bunch of parameters defining the properties of the synchrotron power spectra, and of synchrotron Spectral Energy Distribution (SED), the output are the stokes IQU maps simulated as Gaussian random fields of the defined spectra.
In particular, synchrotron power spectra :math:`C_{\ell}` are assumed to follow a power law as a function of :math:`\ell`: :math:`C_{\ell}^{TT/TE/EE/BB}\propto\ell^{\alpha}`.
Spectra are defined by:

1. The slope :math:`\alpha` (same for all the spectra)
2. The amplitude of TT and EE spectra at :math:`\ell=80`,
3. The ratio between B and E-modes

Stokes Q and U maps are generated as random realization of the polarization spectra. For the temperature map the situation is slightly different as we want the total intensity map to be positive everywhere.
The Stokes I map is generated in the following way:

if target Nside<=64:
    1. The TT power spectrum is  :math:`C_\ell \propto \ell^\alpha` and :math:`C_\ell[0]=0`
    2. A first temparature map T is generated as a gaussian realization of this power spectrum
    3. A new map is obtained by adding to T an offset whose value is taken from a reference map
    4. If T+offset is positive everywhere than this is the output temperature map
    5. Otherwise a cut in the TT power spectrum is applied in the following way: :math:`C_\ell[1:\ell_{cut}] = C_\ell[\ell_{cut}]`
    6. A new T+offset map is generated. The value of :math:`\ell_{cut}` is the minimum one for which T+offset is positive everywhere

if target Nside>64:
    1. a map at Nside=64 is generated following the procedure above and then filtered to retain only large angular scales (ell<30)
    2. a map at the target Nside is generated including only small scales (ell>30) with the same seed as the map at point 1.
    3. the two maps are added together
    4. In case the coadded map still has negative pixels a small offset is added to make it positive everywhere


The default parameters are optimized for SO-SAT observations. Meaning that the amplitudes of power spectra are normalized in the 10% sky region observed by the instrument. In particular:

1. The amplitude of TT spectrum is taken from PySM-s0 model at 23GHz.

   TT_amplitude = 20 :math:`\mu K^2` (for :math:`D_\ell` at :math:`\ell=80`)
2. The offset for T map is also taken from PySM-s0 model at 23GHz.

   Toffset = 72 :math:`\mu K`
2. The amplitude of EE spectrum is taken from S-PASS at 2.3GHz extrapolated at 23GHz with a powerlaw with :math:`\beta_s=-3.1`

   EE_amplitude = 4.3 math:`\mu K^2` (for :math:`D_\ell` at :math:`\ell=80`)
3. ratio between B and E modes from Krachmalnicoff et al. 2018

   B_to_E = 0.5
4. spectral tilt from Krachmalnicoff et al 2018

   alpha = -1
5. spectral index from Planck IX 2018

   beta = -3.1
6. Default value for curvature is zero


GaussianDust
===========

This class implements Gaussian simulations for Galactic thermal dust emission.
The inputs are a bunch of parameters defining the properties of dust power spectra, and of dust Spectral Energy Distribution (SED), the output are the stokes IQU maps simulated as Gaussian random fields of the defined spectra.
In particular, dust power spectra :math:`C_{\ell}` are assumed to follow a power law as a function of :math:`\ell`: :math:`C_{\ell}^{TT/TE/EE/BB}\propto\ell^{\alpha}`.
Spectra are defined by:

1. The slope :math:`\alpha` (same for all the spectra)
2. The amplitude of TT and EE spectra at :math:`\ell=80`,
3. The ratio between B and E-modes
4. The degree of correlation between T and E.

Stokes Q and U maps are generated as random realization of the polarization spectra. For the temperature map the situation is slightly different as we want the total intensity map to be positive everywhere.
The Stokes I map is generated in the following way:

if target Nside<=64:
    1. The TT power spectrum is  :math:`C_\ell \propto \ell^\alpha` and :math:`C_\ell[0]=0`
    2. A first temparature map T is generated as a gaussian realization of this power spectrum
    3. A new map is obtained by adding to T an offset whose value is taken from a reference map
    4. If T+offset is positive everywhere than this is the output temperature map
    5. Otherwise a cut in the TT power spectrum is applied in the following way: :math:`C_\ell[1:\ell_{cut}] = C_\ell[\ell_{cut}]`
    6. A new T+offset map is generated. The value of :math:`\ell_{cut}` is the minimum one for which T+offset is positive everywhere

if target Nside>64:
    1. a map at Nside=64 is generated following the procedure above and then filtered to retain only large angular scales (ell<30)
    2. a map at the target Nside is generated including only small scales (ell>30) with the same seed as the map at point 1.
    3. the two maps are added together
    4. In case the coadded map still has negative pixels a small offset is added to make it positive everywhere

Typical values for :math:`\ell_{cut}` are between :math:`\ell=4` and :math:`\ell=9`, depending on realization (and also on the Nside of the output map). This implementation removes some power at the very large scales.

The default parameters are optimized for SO-SAT observations. Meaning that the amplitudes of power spectra are normalized in the 10% sky region observed by the instrument. In particular:

1. The amplitude of TT spectrum is taken from PySM-d0 model at 353GHz.

   TT_amplitude = 350 :math:`\mu K^2` (for :math:`D_\ell` at :math:`\ell=80`)
2. The offset for T map is also taken from PySM-d0 model at 353GHz.

   Toffset = 18 :math:`\mu K`
2. The amplitude of EE spectrum is taken from Planck map at 353GHz

   EE_amplitude = 100 math:`\mu K^2` (for :math:`D_\ell` at :math:`\ell=80`)
3. ratio between B and E modes from Planck IX 2018

   B_to_E = 0.5
4. spectral tilt from Planck IX 2018

   alpha = -0.42
5. spectral index and temperature from Planck IX 2018

   beta = 1.53, T=19.6 K
