.. _cog:

###############
Curve of Growth
###############

The curve of growth method (CoG) is described in `Yu et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...898..102Y/abstract>`_.
``dysh`` tries to follow the description given by `Yu et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...898..102Y/abstract>`_, but there are some differences.
Here we describe the ``dysh`` implementation of the CoG method, which is available through the :py:meth:`Spectrum.cog <dysh.spectra.spectrum.Spectrum.cog>` method or through the :py:func:`dysh.spectra.core.curve_of_growth`.
An example of how to use :py:meth:`Spectrum.cog <dysh.spectra.spectrum.Spectrum.cog>` is provided in the :doc:`HI survey tutorial <users_guide/hi_survey.html>`.

Central Velocity
================

If no central velocity is specified (parameter ``vc``), it will be computed using the first moment of the spectrum inside the ranges defined by ``bchan`` and ``echan``.
That is

.. math::

    v_{\mathrm{c}}=\frac{\sum_{i} T(v_{i})v_{i}}{\sum_{i} T(v_{i})},

where :math:`T(v_{i})` are the flux values at channels :math:`v_{i}`.

This method is not robust, meaning that small baseline deviations from a flat response will bias the result.
To minimize these effects we recommend selecting a spectral range around the spectral line of interest, either using the ``bchan`` and ``echan`` arguments or by cropping the spectrum.

If no value of ``vc`` is provided, the uncertainty in the estimated :math:`v_{\mathrm{c}}` is computed as

.. math::

    \sigma_{v_{\mathrm{c}}}^{2}=\left(\frac{v_{\mathrm{c}}}{\sum_{i} T(v_{i})v_{i}}\sqrt{\sum (v\sigma)^{2}}\right)^{2}+\left(\sqrt{N}\sigma \frac{v_{\mathrm{c}}}{\sum_{i} T(v_{i})}\right)^{2},

where :math:`N` is the number of channels and :math:`\sigma` is the rms in the line free channels.

Line Area
=========

For the line area (or line intensity, flux intensity, integrated intensity) we use the same definition as `Yu et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...898..102Y/abstract>`_, their equation (1).
To determine the line area, :math:`F_{\mathrm{t}}`, we take the median of :math:`F_{\mathrm{t}}(v)` after it becomes flat.
To determine the point at which :math:`F_{\mathrm{t}}(v)` becomes flat, we estimate the slope of :math:`F_{\mathrm{t}}(v)` and take the point at which the slope, :math:`s`, satisfies

.. math::

    s<f\sigma_{s},

with :math:`f` being the parameter ``flat_tol`` (defaults to 0.1) and :math:`\sigma_{s}` is the rms of :math:`s`.
The estimated line area is returned as the ``flux`` entry in the return dictionary.

The uncertainty in the line area is estimated as the rms of :math:`F_{\mathrm{t}}(v)` after it becomes flat.
This also incorporates 3% of the line area, which was determined empirically using synthetic spectra (added in quadrature to the rms of :math:`F_{\mathrm{t}}(v)`).
The uncertainty in the line area is returned as the ``flux_std`` entry in the return dictionary.

Line Width
==========

For the line width we use the same definition as `Yu et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...898..102Y/abstract>`_.
As the authors note, this is different than :math:`W_{x}` in that the line widths are those that enclose a certain fraction of the line intensity, whereas :math:`W_{x}` is the line width at a fraction :math:`x` of the peak flux.
The line width is returned as the ``width`` entry in the return dictionary.
This entry is another dictionary, where each key is the fraction of the flux and the value the width at that fraction.
For example,

.. code:: python

    cog = Spectrum.cog()
    cog["width"][0.5]

is the line width that encompases 50% of the line area.

To estimate the uncertainty in the line width we compute the error in the normalized curve of growth, :math:`\hat{F}=F_{\mathrm{t}}(v)/F_{\mathrm{t}}`,

.. math::

    \sigma_{\hat{F}}^{2}=\left(\frac{\hat{F}}{F_{\mathrm{t}}(v)}\sigma_{F_{\mathrm{t}}(v)}\right)^{2}+\left(\frac{\hat{F}}{F_{\mathrm{t}}}\sigma_{F_{\mathrm{t}}}\right)^{2}

with :math:`\sigma_{F_{\mathrm{t}}(v)}` the uncertainty in the curve of growth :math:`F_{\mathrm{t}}(v)` and :math:`\sigma_{F_{\mathrm{t}}}` the uncertainty in the line area.
We estimate :math:`\sigma_{F_{\mathrm{t}}(v)}` as the rms in the line free channels times the channel width.
Then, we compute the width adding and subtracting :math:`\sigma_{\hat{F}}` to :math:`\hat{F}` at each fraction of the line area.
The final uncertainty is the maximum between the difference of the width and the width plus :math:`\sigma_{\hat{F}}`, the width minus :math:`\sigma_{\hat{F}}`, and the channel width.
We add 1% of the line width to the uncertainty in the line width, which was determined empirically using synthetic spectra.
The uncertainty in the line width is returned as the ``width_std`` entry in the return dictionary.
And, as the ``width`` entry, it is also a dictionary.

Flux and Shape Assymetry
========================

For the flux assymetry, :math:`A_{F}`, we use the definition of `Yu et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...898..102Y/abstract>`_.
This is returned as the ``A_F`` entry in the return dictionary.
We do not provide an estimate for the uncertainty in :math:`A_{F}`.

For the second assymetry parameter, :math:`A_{C}`, we use the definition of `Yu et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...898..102Y/abstract>`_.
This is returned as the ``A_C`` entry in the return dictionary.
We do not provide an estimate for the uncertainty in :math:`A_{C}`.

The assymetry parameters :math:`A_{F}` and :math:`A_{C}` are close to unity for symmetric line profiles, and greater for more assymetric line profiles.

Concentration
=============

For the concentration of the line profile, :math:`C_{V}`, we use the definition of `Yu et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...898..102Y/abstract>`_

.. math::

    C_{V}=\frac{V_{85}}{V_{25}}

where :math:`V_{85}` and :math:`V_{25}` are the line widths at 85% and 25% of the total flux, respectively.
The concentration is returned as the ``C_V`` parameter in the return dictionary.
We do not provide an estimate for the uncertainty in :math:`C_{V}`.

As noted by `Yu et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...898..102Y/abstract>`_, a Gaussian profile has :math:`C_{V}=3.9`, while a boxcar prifile has :math:`C_{V}=3.4`.
