.. _sdmath:


Single Dish Math
~~~~~~~~~~~~~~~~

Here we briefly review some of the equations governing Single Dish math.


Getting the system temperature


.. math::  :label: eq_sdmath1

   T_{sys} = T_{amb} { { SKY } \over { HOT - SKY } }

or

.. math:: :label: eq_sdmath2

   T_{sys} = T_{cal} { { <SKY> } \over { <HOT - SKY> } } + T_{cal}/2

where the :math:`< >` operator averages over the center (typically 80%) portion of the spectrum.
This way :math:`T_{sys}` is a scalar. The routine ``meantsys`` computes this.  The HOT and SKY
are also referred to sometimes as CAL and SIG.

.. math:: :label: eq_sdmath3

   T_A = T_{sys}  {   { SIG - REF } \over {REF} }

All of these have values for each channel. How exactly the :math:`T_{sys}` is computed (scalar, vector,
mean/median) can vary with different implementations.

Note in some places you may see the SIG/REF referred to as ON/OFF.	  



The radiometer equation equates the noise to the system temperature, integration time and bandwidth

.. math:: :label: eq_radiometer


   \Delta T =  \frac{T_{\rm{sys}}}{\sqrt{\Delta t \Delta\nu}}


The effective beam (see also GBT memo 296, and gbtpipe/Gridding.py - at 109 GHz)

..  1.18 * (c / nu0 / 100.0) * 180 / np.pi  # in degrees

.. math:: :label: eq_beam109

    1.18 {\lambda \over D}

with D=100 m.


Reduction of noise with **smoothref=N**:

.. math:: :label: eq_smoothref

     \sigma_N = \sigma_1 \sqrt{   {N+1} \over  {2N}  }

Weight factors are 1/:math:`\sqrt(RMS)` following the radiometer equation

.. math:: :label: eq_sdmath5

      \frac{\Delta t \Delta\nu}{T_{\rm{sys}}^{2}}

Effective exposure time in an ON/OFF observation

.. math:: :label: eq_sdmath6

  { t_{ON} * t_{OFF} } \over {  t_{ON} + t_{OFF}  }

The Ruze equation equates
the gain of an antenna to the root mean square (:math:`\epsilon`) of the antenna's random surface errors.

.. math:: :label: eq_ruze

   G = G_0 \exp{ (-4\pi\epsilon / \lambda^2) }


As shown in :eq:`eq_sdmath2` we can ...


Something about Doppler and Velocity Frames?


See also  :ref:`cog` for math behind the Curve of Growth method.


Temperature scales:   Ta, Ta', Ta*, Tmb -
Correcting for Atmospheric Opacity -
see https://library.nrao.edu/public/memos/gbt/GBT_302.pdf
