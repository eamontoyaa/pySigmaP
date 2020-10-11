.. highlight:: shell

============
Installation
============


Stable release
--------------

To install **pySigmaP**, run this command in your terminal:

.. code-block:: console

    $ pip install pysigmap

This is the preferred method to install pySigmaP, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for pySigmaP can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/eamontoyaa/pysigmap

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/eamontoyaa/pysigmap
.. _tarball: https://github.com/eamontoyaa/pysigmap/tarball/master


Requirements
------------

The code was written in Python 3. The packages `numpy <http://www.numpy.org/>`_,
`scipy <https://www.scipy.org/>`_, `matplotlib <https://matplotlib.org/>`_, 
`pandas <https://pandas.pydata.org/>`_
and `scikit-learn <https://scikit-learn.org/stable/>`_ are
required for using **pySigmaP**. All of them are
downloadable from the PyPI repository by opening a terminal and typing the
following code lines:


::

    pip install numpy
    pip install scipy
    pip install matplotlib
    pip install pandas
    pip install scikit-learn

`LaTeX <https://www.latex-project.org/>`_ is also required to render fonts. Linux distributions
usually have pre-installed LaTeX; however, for Microsoft Windows operative systems, it is suggested
to install `MiKTeX <https://miktex.org/download/>`_,  an up-to-date implementation of TeX/LaTeX.

