"""
pachecosilva.py module.

Contains the class and its methods for determinig the preconsolidation
pressure from a consolidation test by the method proposed by
Pacheco Silva (1970).

References
----------
Pacheco Silva, F. 1970. A new graphical construction for determination of the
pre-consolidation stress of a soil sample. In Proceedings of the 4th Brazilian
Conference on Soil Mechanics and Foundation Engineering, Rio de Janeiro,
Brazil. Vol. 2, No.1, pp. 225–232.

"""

# -- Required modules
import numpy as np
from numpy.polynomial.polynomial import polyfit, polyval
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import matplotlib.ticker as mtick

plt.rcParams['font.family'] = 'Serif'
plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True


class PachecoSilva():
    """PachecoSilva class."""

    def __init__(self, data):
        """
        Initialize the PachecoSilva class.

        Instance an object to perform the PachecoSilva's method for determining
        the preconsolidation pressure from a unidimensional consolidation test.

        Parameters
        ----------
        data : Object instanced from the Data class.
            Contains the data structure from the consolidation test.

        Returns
        -------
        None.

        """
        self.data = data
        return

    def getSigmaP(self, range2fitTOP=None, range2fitNCL=None):
        """
        Return the value of the preconsolidation pressure or yield stress

        Parameters
        ----------
        range2fitTOP : list, tuple or array (length=2), optional
            Initial and final pressures between which the third-order
            polynomial (TOP) is fit to the compressibility curve. If None, the
            TOP is fit to the second third of the curve. The default is None.
        range2fitNCL : list, tuple or array (length=2), optional
            Initial and final pressures between which the first-order
            polynomial is fit to the compressibility curve on the normally
            consolidated line (NCL). If None, the NCL is fit to the last three
            points of the curve. The default is None.

        Returns
        -------
        fig : matplotlib figure
            Figure with the development of the method and the results.

        """
        if range2fitTOP is None:  # Indices for fitting the third order poly
            idxInitTOP = int(np.floor(len(self.data.cleaned) / 3))
            idxEndTOP = int(np.ceil(2 * len(self.data.cleaned) / 3))
        else:
            idxInitTOP = self.data.findStressIdx(
                stress2find=range2fitTOP[0], cleanedData=True)
            idxEndTOP = self.data.findStressIdx(
                stress2find=range2fitTOP[1], cleanedData=True)

        if range2fitNCL is None:  # Indices for fitting the NCL line
            idxInitNCL, idxEndNCL = -3, None
        else:
            idxInitNCL = self.data.findStressIdx(
                stress2find=range2fitNCL[0], cleanedData=True)
            idxEndNCL = self.data.findStressIdx(
                stress2find=range2fitNCL[1], cleanedData=True) if \
                range2fitNCL[1] < self.data.cleaned['stress'].max() else None

        # -- fittig a third order polynomial to data without unloads
        sigmaTOP = self.data.cleaned['stress'][idxInitTOP: idxEndTOP]
        sigmaTOPlog = np.log10(sigmaTOP)
        eTOP = self.data.cleaned['e'][idxInitTOP: idxEndTOP]
        p1_0, p1_1, p1_2, p1_3 = polyfit(sigmaTOPlog, eTOP, deg=3)
        r2TOP = r2_score(
            y_true=eTOP, y_pred=polyval(sigmaTOPlog, [p1_0, p1_1, p1_2, p1_3]))
        xFitTOP = np.linspace(sigmaTOP.iloc[0], sigmaTOP.iloc[-1], 500)
        yFitTOP = polyval(np.log10(xFitTOP), [p1_0, p1_1, p1_2, p1_3])

        # -- Linear regresion of points on normally consolidated line (NCL)
        sigmaNCL = self.data.cleaned['stress'][idxInitNCL: idxEndNCL]
        sigmaNCLlog = np.log10(sigmaNCL)
        eNCL = self.data.cleaned['e'][idxInitNCL: idxEndNCL]
        p2_0, p2_1 = polyfit(sigmaNCLlog, eNCL, deg=1)
        r2NCL = r2_score(
            y_true=eNCL, y_pred=polyval(sigmaNCLlog, [p2_0, p2_1]))

        # -- Lines of the Pacheco Silva's method
        # Intersection NCL - e_0
        xPt1 = 10 ** ((self.data.e_0 - p2_0) / p2_1)
        xFitNCL = np.linspace(xPt1, sigmaNCL.iloc[-1])
        yFitNCL = polyval(np.log10(xFitNCL), [p2_0, p2_1])
        # Intersection first vertical - compressibility curve
        yPt2 = polyval(np.log10(xPt1), [p1_0, p1_1, p1_2, p1_3])
        # Intersection first horizontal - NCL (Preconsolidation pressure)
        self.sigmaP = 10 ** ((yPt2 - p2_0) / p2_1)

        # -- plotting
        fig = plt.figure(figsize=[9, 4.8])
        # ax = fig.add_subplot(111)  # Compressibility curve
        ax = fig.add_axes([0.08, 0.12, 0.55, 0.85])
        ax.plot(self.data.raw['stress'][1:], self.data.raw['e'][1:], ls='--',
                marker='o', lw=1, c='k', mfc='w',
                label='Compressibility curve')
        ax.plot(sigmaNCL, eNCL, ls='', marker='.', lw=0.8, color='crimson')
        ax.plot(xFitNCL, yFitNCL, ls='--', lw=0.8, color='crimson',
                label=f'NCL linear fit (R$^2={r2NCL:.3f}$)')
        ax.plot(sigmaTOP, eTOP, ls='', marker='.', lw=0.8, color='darkcyan')
        ax.plot(xFitTOP, yFitTOP, ls='--', lw=0.8, color='darkcyan',
                label=str().join([r'$3^\mathrm{rd}$-order polynomial fit',
                                  f'\n(R$^2={r2TOP:.3f}$)']))  # poly fit
        # Lines of the Pacheco-Silva's method
        ax.hlines(y=self.data.e_0, xmin=0, xmax=xPt1, lw=0.8, ls='-',
                  color='k', label=f'$e_0 = {self.data.e_0:.2f}$')
        ax.vlines(x=xPt1, ymin=yPt2, ymax=self.data.e_0, lw=0.8, color='k',
                  ls='-.')
        ax.hlines(y=yPt2, xmin=xPt1, xmax=self.sigmaP, lw=0.8, color='k',
                  ls='-.')
        ax.vlines(x=self.sigmaP, ymin=yFitNCL.min(), ymax=yPt2, lw=0.8,
                  color='k', ls='-.')
        ax.plot(self.sigmaP, yPt2, ls='', marker='D', c='r', ms=5, mfc='w',
                label=str().join([r'$\sigma^\prime_\mathrm{p}=$ ',
                                  f'{self.sigmaP:.0f} kPa']))
        # other details
        ax.set(xscale='log', ylabel='Void ratio $(e)$',
               xlabel=str().join(['Vertical effective stress ',
                                  r'$(\sigma^\prime_\mathrm{v})$ [kPa]']))
        ax.xaxis.set_major_formatter(mtick.ScalarFormatter())
        ax.grid(True, which="both", ls='--', lw=0.5)
        ax.legend(bbox_to_anchor=(1.04, 0.5), loc=6,
                  title=r"\textbf{Pacheco-Silva's method}")
        return fig


# %%
"""
BSD 2 license.

Copyright (c) 2020, Exneyder A. Montoya-Araque and Alan J. Aparicio-Ortube.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
