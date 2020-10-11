========
pySigmaP
========

.. image:: https://img.shields.io/badge/Made%20with-Python3-brightgreen.svg
        :target: https://www.python.org/
        :alt: made-with-python

.. image:: https://img.shields.io/pypi/v/pysigmap.svg
        :target: https://pypi.python.org/pypi/pysigmap

.. image:: https://img.shields.io/badge/License-BSD%202--Clause-brightgreen.svg
        :target: https://github.com/eamontoyaa/pysigmap/blob/master/LICENSE
        :alt: License

.. image:: https://readthedocs.org/projects/pysigmap/badge/?version=latest
        :target: https://pysigmap.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Open-source application software developed in Python3 for interpreting
the preconsolidation pressure of fine-grained soils in incremental loading
oedometer testing. pySigmaP includes nine methods such as the methods of Casagrande,
Pacheco-Silva, Butterfield, Oikawa, Becker et al., Morin, Onitsuka et al.,
Wang and Frost, and Boone.


* Free software: BSD 2-Clause License
* Documentation: https://pysigmap.readthedocs.io.


Features
--------

* `Documentation <https://pysigmap.readthedocs.io>`_
* `PyPI <https://pypi.org/project/pysigmap>`_
* `GitHub <https://github.com/eamontoyaa/pysigmap>`_
* Open source and free software: `BSD-2-Clause <https://opensource.org/licenses/BSD-2-Clause>`_.


Structure
---------

**pyBIMstab** was written under the object-oriented paradigm and is divided into six modules,
one of them contains the class to load and manage the data of the compressibility curve, the other
five contain the classes that solve the nine methods for getting the preconsolidation pressure
or yield stress.

Each module is well documented and contains basic usage examples.

.. toctree::
   :maxdepth: 10

   data
   casagrande
   energy
   bilog
   pachecosilva
   boone

   


