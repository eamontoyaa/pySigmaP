"""
``boone.py`` module.

Contains the class and its methods for interpreting the preconsolidation
pressure from a compressibility curve via the method proposed by Boone (2010).

References
----------
Boone, S. J. (2010). A critical reappraisal of "preconsolidation
pressure" interpretations using the oedometer test. Canadian Geotechnical
Journal, 47, 3, 281-296. https://doi.org/10.1139/T09-093.

"""

# -- Required modules
import numpy as np
from numpy.polynomial.polynomial import polyval
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

plt.rcParams['font.family'] = 'Serif'
plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True
# High-contrast qualitative colour scheme
colors = ('#DDAA33',  # yellow
          '#BB5566',  # red
          '#004488')  # blue


class Boone():
    """
    ``Boone`` class.

    When the object is instanced, the method ``getSigmaP()`` calculates the
    preconsolidation pressure by the method proposed by Boone (2010) based on
    the compression and recompression indices obtained with the ``Data`` class.

    Attributes
    ----------
    data : Object instanced from the ``Data`` class.
        Contains the data structure from the oedometer test. See the class
        documentation for more information.

    Examples
    --------
    >>> urlCSV = ''.join(['https://raw.githubusercontent.com/eamontoyaa/',
    >>>                   'data4testing/main/pysigmap/testData.csv'])
    >>> data = Data(pd.read_csv(urlCSV), sigmaV=75)
    >>> method = Boone(data)
    >>> method.getSigmaP()
    >>> method.sigmaP, method.ocr
    (384.38457280143547, 5.125127637352473)
    """

    def __init__(self, data):
        """Initialize the class."""
        self.data = data
        return

    def getSigmaP(self):
        """
        Return the value of the preconsolidation pressure.

        Returns
        -------
        fig : matplotlib figure
            Figure with the development of the method and the results.
        """
        # Intercept of the line parallel to Cr passing through (σ_v0, e_σ_v0)
        idxCr2Int = self.data.eSigmaV + \
            self.data.idxCr * np.log10(self.data.sigmaV)
        # Intersection of Line parallel to Cr2 - Cc (Preconsolidation pressure)
        self.sigmaP = 10 ** ((self.data.idxCcInt - idxCr2Int) /
                             (-self.data.idxCr + self.data.idxCc))
        self.eSigmaP = polyval(np.log10(self.sigmaP),
                               [self.data.idxCcInt, -self.data.idxCc])
        self.ocr = self.sigmaP / self.data.sigmaV

        # -- plotting
        fig = plt.figure(figsize=[9, 4.8])
        ax = fig.add_axes([0.08, 0.12, 0.55, 0.85])
        ax.plot(self.data.raw['stress'][1:], self.data.raw['e'][1:],
                ls=(0, (1, 1)), marker='o', lw=1.5, c='k', mfc='w',
                label='Compressibility curve')
        # Compression index (Cc)
        x4Cc = np.linspace(self.sigmaP, self.data.cleaned['stress'].iloc[-1])
        y4Cc = -self.data.idxCc * np.log10(x4Cc) + self.data.idxCcInt
        ax.plot(x4Cc, y4Cc, ls='-', lw=1.125, color=colors[1],
                label=str().join([r'$C_\mathrm{c}=$',
                                  f'{self.data.idxCc:.3f}']))
        if self.data.fitCc:
            ax.plot(self.data.cleaned['stress'].iloc[self.data.maskCc],
                    self.data.cleaned['e'].iloc[self.data.maskCc],
                    ls='', marker='x', color=colors[1],
                    label=f'Data for linear fit\n(R$^2={self.data.r2Cc:.3f}$)')
        # Recompression index (Cr)
        x4Cr = np.linspace(
            self.data.raw['stress'].iloc[self.data.maskCr].min(),
            self.data.raw['stress'].iloc[self.data.maskCr].max())
        y4Cr = -self.data.idxCr * np.log10(x4Cr) + self.data.idxCrInt
        ax.plot(x4Cr, y4Cr, ls='-', lw=1.125, color=colors[2],
                label=str().join([r'$C_\mathrm{r}=$',
                                  f'{self.data.idxCr:.3f}']))
        ax.plot(self.data.raw['stress'].iloc[self.data.maskCr],
                self.data.raw['e'].iloc[self.data.maskCr], ls='', marker='+',
                color=colors[2],
                label=f'Data for linear fit\n(R$^2={self.data.r2Cr:.3f}$)')
        # Line pararel to Cr
        x4Cr2 = np.linspace(
            self.data.cleaned['stress'].iloc[1], self.sigmaP)
        y4Cr2 = polyval(np.log10(x4Cr2), [idxCr2Int, -self.data.idxCr])
        ax.plot(x4Cr2, y4Cr2, ls='--', c=colors[2], lw=1.125,
                label='Parallel line to $C_\\mathrm{r}$')
        # Other plots
        ax.plot(self.data.sigmaV, self.data.eSigmaV, ls='', marker='|', c='r',
                ms=15, mfc='w', mew=1.5,
                label=str().join([r'$\sigma^\prime_\mathrm{v_0}=$ ',
                                  f'{self.data.sigmaV:.0f} kPa']))
        ax.plot(self.sigmaP, self.eSigmaP, ls='', marker='o', c=colors[0],
                ms=7, mfc='w', mew=1.5,
                label=str().join([r'$\sigma^\prime_\mathrm{p}=$ ',
                                  f'{self.sigmaP:.0f} kPa\n',
                                  f'OCR= {self.ocr:.1f}']))
        # Other details
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set(xscale='log', ylabel='Void ratio, $e$',
               xlabel=str().join(['Vertical effective stress, ',
                                  r'$\sigma^\prime_\mathrm{v}$ [kPa]']))
        ax.xaxis.set_major_formatter(mtick.ScalarFormatter())
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.grid(False)
        ax.legend(bbox_to_anchor=(1.125, 0.5), loc=6,
                  title=r"\textbf{Boone method}")
        return fig


# %%
"""
2-Clause BSD License.

Copyright 2020, EXNEYDER A. MONTOYA-ARAQUE, A. J. APARICIO-ORTUBE,
DAVID G. ZAPATA-MEDINA, LUIS G. ARBOLEDA-MONSALVE AND
UNIVERSIDAD NACIONAL DE COLOMBIA.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
