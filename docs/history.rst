=======
History
=======

0.1.0 (2020-10-10)
------------------

* First release on PyPI.

0.1.1 (2020-10-10)
------------------

* Minor updates to html documentation.

0.1.2 (2020-10-11)
------------------

* Minor updates to html documentation.

0.1.3 (2020-11-30)
------------------

* Improvements to figures styles.

0.1.4 (2020-11-21)
------------------

* Minor updates to html documentation.

0.1.5 (2020-12-26)
------------------

* Minor updates to html documentation.

0.1.6 (2020-12-26)
------------------

* Minor updates to html documentation and improvements to figures.

0.1.7 (2020-12-26)
------------------

* Use mstools instead of scikit-learn for linear regression coefficient of determination.

0.1.8 (2021-06-28)
------------------

* Minor improvements to figures.

0.1.9 (2022-08-24)
------------------

* Update documentation.
* Include feature to the ``Data`` class to read data with multiple unloading-reloading stages.
* Remove dependency of external LaTeX installation.

0.1.10 (2023-10-31)
-------------------

* Include the ``range2fitCS`` to the ``getSigmaP`` method of the ``Casagrande`` class to limit the stress range for the cubic spline in the calculation of the maximum curvature point.
* The maximum curvature point is now calculated using a function within the ``casagrande`` module
  that determines the maximum value of the curvature that is not at the ends, instead
  of the ``find_peaks`` SciPy's function. However, if the value is not determined, the maximum
  curvature point is calculated as the absolute maximum value of the curvature.
* Improve some figures:

    - The void ratio of the :math:`\sigma'_\mathrm{v}` value is now determined by
      projecting :math:`\sigma'_\mathrm{v}` over a linear interpolation instead of a cubic spline.
    - :math:`C_\mathrm{c}` and :math:`C_\mathrm{r}` are not calculated inmediatelly when loading the data.
      Instead, they are calculated when theis respective methods are called.
    - Other minor improvements.
