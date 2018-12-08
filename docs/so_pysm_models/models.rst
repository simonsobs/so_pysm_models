Summary of Models
**********************

Synchrotron
===========

The difference between stokes Q/U and I is that we want the total intensity map to be positive everywhere, meaning that the same implementation used for Q/U cannot be applied in this case.

The way I bypass this problem is the following:

1. I generate a power-law TT power spectrum: :math:`C_\ell \propto \ell^\alpha`
2. I put :math:`C_\ell[0]=0`
3. I generate a temperature map T as a gaussian realization of this power spectrum
4. I add at this T map an offset whose value is taken from a reference map (in this case the offset is the mean value of the pysm T map at 23 GHz in the SO-SAT region)
5. I check whether the T+offset map is positive everywhere
6. If not I applied a cut in the TT power spectrum in the following way: :math:`C_\ell[1:\ell_{cut}] = C_\ell[\ell_{cut}]` and generate again the T+offset map. The value of lcut is the minimum one for which T+offset is positive everywhere

Typical values for lcut are between :math:`\ell=4` and :math:`\ell=9`, depending on realization (and also on the Nside of the output map).
This implementation removes some power at the very large scales, which nevertheless are nor probed by SO.

The code is slower now, as it has to generate the output map several times, until T is positive everywhere.
