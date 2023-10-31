# `pySigmaP`

[![made-with-python](https://img.shields.io/badge/Made_with-Python3-306998?style=badge&logo=python&logoColor=white)](https://github.com/eamontoyaa/pySigmaP)      [![PyPI repo](https://img.shields.io/pypi/v/pysigmap.svg)](https://pypi.org/project/pysigmap)      [![BSD-2 License](https://img.shields.io/badge/License-BSD2-FFAB00?style=badge&logo=opensourceinitiative&logoColor=white)](https://opensource.org/license/bsd-2-clause/)            [![Docs](https://readthedocs.org/projects/pysigmap/badge/?version=latest)](https://pysigmap.readthedocs.io/en/latest/?badge=latest)

## Description

**``pySigmaP``** is an application software developed in Python3 to determine
the preconsolidation pressure of soils in incremental loading (IL) oedometer
testing. **``pySigmaP``** includes nine methods such as the method of Casagrande,
Pacheco Silva, Butterfield, Oikawa, Becker et al., Morin, Onitsuka et al.,
Wang and Frost, and Boone.

In this repo and in the package [documentation](https://pysigmap.readthedocs.io) you will find the source code, instructions for installation, docstrings, examples of use, development history track, and references.

## Installation

It is suggested to create a virtual environment to install and use the program.

### Stable release

We recommend installing **``pySigmaP``** from [PyPI](https://pypi.org/project/pySigmaP), as it will always install the most recent stable release.  To do so, run this command in your terminal:

    pip install pysigmap

### From sources

The sources for **``pySigmaP``** can be downloaded from the [Github repo](https://github.com/eamontoyaa/pySigmaP). You can clone the public repository running the following command:

    git clone git://github.com/eamontoyaa/pySigmaP

Once you have a copy of the source code, you can install it with the following instruction:

    pip install -e .

### Dependencies

The code was written in Python 3.9. The packages `numpy`, `matplotlib`, `scipy`, `pandas` and `mstools` are required for using **``pySigmaP``**. They should be installed along with the package, however, all of them can also be manually installed from the PyPI repository by opening a terminal and typing the following code lines:

    pip install numpy==1.23.5
    pip install matplotlib==3.7.1
    pip install scipy==1.11.3
    pip install pandas==1.5.3
    pip install mstools==0.1.0


## Authorship and Citation

The team behind **``pySigmaP``** includes:

**Development Team**

* Exneyder A. Montoya-Araque, Geol.Eng., MEng., <eamontoyaa@unal.edu.co>
* Alan J. Aparicio-Ortube, Civ.Eng., MSc., PhD. Student,  <aaparicioo@unal.edu.co>

**Contributors/Advisors**

* David G. Zapata-Medina, Civ.Eng., MSc., PhD. <dgzapata@unal.edu.co>
* Luis G. Arboleda-Monsalve, Civ.Eng., MSc., PhD. <luis.arboleda@ucf.edu>


To cite **``pySigmaP``** in publications, use the following reference:

    Montoya-Araque et al. (2022). An open-source application software to determine the preconsolidation pressure of soils in incremental loading oedometer testing: pySigmaP. SoftwareX, 17, 100990. https://doi.org/10.1016/j.softx.2022.100990

A BibTeX entry for LaTeX users is:

``` bibtex
@article{MontoyaAraque_etal_2022_art,
    author = {Montoya-Araque, Exneyder A. and Aparicio-Ortube, A.J. and Zapata-Medina, David G. and Arboleda-Monsalve, Luis G.},
    doi = {10.1016/j.softx.2022.100990},
    issn = {23527110},
    journal = {SoftwareX},
    keywords = {Oedometer testing,Open-source application software,Preconsolidation pressure,Python 3,Soil},
    month = {jan},
    pages = {100990},
    publisher = {Elsevier B.V.},
    title = {{An open-source application software to determine the preconsolidation pressure of soils in incremental loading oedometer testing: pySigmaP}},
    volume = {17},
    year = {2022}
}
```

## License

**``pySigmaP``** is a free and open-source application sorfware licensed under the terms of the [BSD-2-Clause License](https://opensource.org/licenses/BSD-2-Clause).
