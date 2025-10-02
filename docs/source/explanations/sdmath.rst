.. _sdmath:


Single Dish Math
~~~~~~~~~~~~~~~~

The meat of Single Dish math is getting the system temperature


.. math::  :label: sdmath1

   T_{sys} = T_{amb} { { SKY } \over { HOT - SKY } }

or

.. math:: :label: sdmath2

   T_{sys} = T_{cal} { { <SKY> } \over { <HOT - SKY> } } + T_{cal}/2

where the :math:`< >` operator averages over the center 80% of the spectrum.
This way :math:`T_{sys}` is a scalar. The routine ``meantsys`` computes this.

and using this system temperature, calculating the signal by comparing an *ON* and *OFF* position,
assuming there is only sky in the *OFF*:

.. math:: :label: sdmath3

   T_A = T_{sys}  {   { ON - OFF } \over {OFF} }

All of these have values for each channel. How exactly the :math:`T_{sys}` is computed (scalar, vector,
mean/median) is something we generally leave open.

.. math::


the effective beam (see also GBT memo 296, and gbtpipe/Gridding.py - at 109 GHz)

.. code-block::

    1.18 * (c / nu0 / 100.0) * 180 / np.pi  # in degrees

Reduction of noise with smoothref=N:

.. math:: :label: sdmath4

     \sigma_N = \sigma_1 \sqrt{   {N+1} \over  {2N}  }

Weight factors are 1/:math:`\sqrt(RMS)`

.. math:: :label: sdmath4

      \frac{\Delta t \Delta\nu}{T_{\rm{sys}}^{2}}

Effective exposure time

.. math:: :label: sdmath5

  { t_{ON} * t_{OFF} } \over {  t_{ON} + t_{OFF}  }

As shown in :eq:`sdmath2` we can ...
