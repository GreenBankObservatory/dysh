.. _sdmath:


Single Dish Math
~~~~~~~~~~~~~~~~

Here we briefly review some of the equations governing Single Dish math.

The Model
=========

The power being recorded by a backend can be described by

.. math::

   P=G\left(T_{\mathrm{A}}+T_{\mathrm{sys}}\right)

with :math:`G` the telescope gain, :math:`T_{\mathrm{A}}` the antenna temperature and :math:`T_{\mathrm{sys}}` the system temperature.

During observations, a noise diode or a hot load can be used to determine the telescope gain or system temperature.
When the noise diode is firing, the power becomes

.. math::

   P^{\mathrm{cal}}=G\left(T_{\mathrm{A}}+T_{\mathrm{sys}}+T_{\mathrm{cal}}\right)

with :math:`T_{\mathrm{cal}}` the equivalent temperature of the noise diode or hot load, a quantity that is known beforehand or that must be derived from calibration observations.

System Temperature
==================

For observations with a noise diode or hot load the system temperature can be computed using

.. math::

   T_{\mathrm{sys}}=T_{\rm{cal}}\left[\frac{P_{\rm{ref}}}{P_{\rm{ref}}^{\rm{cal}}-P_{\rm{ref}}}\right]+\frac{T_{\rm{cal}}}{2}

where :math:`P_{\rm{ref}}` is the power measured towards a reference position, a region without signal.
The factor :math:`T_{\rm{cal}}/2` accounts for the noise diode being fired half of the time — this factor is ommited if using a hot load.
To minimize the "noise" in the computation, dysh takes the average of the numerator and the denominator over the inner 80% of the channels.
So, in practice the system temperature is computed as

.. math::

   T_{\mathrm{sys}}=T_{\rm{cal}}\left[\frac{\langle P_{\rm{ref}}\rangle}{\langle P_{\rm{ref}}^{\rm{cal}}-P_{\rm{ref}}\rangle}\right]+\frac{T_{\rm{cal}}}{2}

where the :math:`\langle\rangle` operator denotes an average. Thus, the system temperatures computed are scalars. In dysh the function responsible for this calculation is :py:func:`dysh.spectra.core.mean_tsys`.
For more details on how to compute :math:`T_{\mathrm{cal}}` for a hot load see `GBT memo #302 <https://library.nrao.edu/public/memos/gbt/GBT_302.pdf>`_.

Antenna Temperature
===================

The antenna temperature is computed using

.. math:: :label: eq_sdmath3

   T_{\rm{A}}=T_{\rm{sys}}\frac{P-P_{\rm{ref}}}{P_{\rm{ref}}}

with :math:`P` the power of the signal (e.g., towards the target in a position switched observation).

Radiometer Equation
===================

The radiometer equation equates the noise to the system temperature, integration time and bandwidth

.. math:: :label: eq_radiometer

   \sigma(T)=\frac{T_{\rm{sys}}}{\sqrt{\Delta t \Delta\nu}}

with :math:`\Delta t` the integration time and :math:`\Delta\nu` the bandwidth.
This is the standard deviation of the signal if it was purely thermal noise.

Weights
-------

By default, dysh uses the inverse variance of the thermal noise as weights

.. math::

   w=\frac{\Delta t \Delta\nu}{T_{\rm{sys}}^{2}}.

Brightness Scales
=================

The definitions of the brightness scales used by dysh are in `GBT memo #302 <https://library.nrao.edu/public/memos/gbt/GBT_302.pdf>`_.


Misc
====

Smoothing the Reference
-----------------------

During calibration it is possible to smooth the reference using the ``smoothref`` argument (this is an argument to the calibration routines, e.g., :py:func:`~dysh.fits.GBTFITSLoad.getps`).
For purely thermal noise this would reduce the noise in the calibrated spectrum by a factor

.. math::

   \sqrt{\frac{N+1}{2N}}

where :math:`N` is the width, in channels, of the smoothing kernel.
So, if using ``smoothref=3`` the noise should be reduced by :math:`\sqrt{4/6}`.

Exposure Time
-------------

Effective exposure time after calibrating using a noisy reference power

.. math:: :label: eq_sdmath6

   \Delta t=\frac{t_{\rm{sig}}t_{\rm{ref}}}{t_{\rm{sig}}+t_{\rm{ref}}}.

Ruze Equation
-------------

The Ruze equation relates the gain of an antenna to the root mean
square (:math:`\delta`) of the antenna's random surface errors.

.. math:: :label: eq_ruze

   G = G_0 \exp{ (-(4\pi\delta / \lambda)^2) }

but the associated beam spreading is a different story.

.. Something about Doppler and Velocity Frames?
