"""
bilog.py module.

Contains the class and its methods for determinig the preconsolidation
pressure from a consolidation test by the bilogarithminc methods proposed by
Butterfield (1979), Oikawa (1987) and Onitsuka etal (1995).

References
----------
Butterfield, R. (1979). A natural compression law for soils (an advance on
e –log p'). Géotechnique, 29, 4, 469-480,
https://doi.org/10.1680/geot.1979.29.4.469

Oikawa, H. (1987). Compression Curve of Soft Soils. Soils and Foundations,
27, 3, 99-104, https://doi.org/10.3208/sandf1972.27.3_99

Onitsuka, K., Hong, Z., Hara, Y., & Yoshitake, S. (1995). Interpretation of
Oedometer Test Data for Natural Clays. Soils and Foundations, 35, 3, 61-70,
https://doi.org/10.3208/sandf.35.61
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


class Bilog():
    """Bilog class."""

    def __init__(self, data):
        """
        Initialize the Bilog class.

        Instance an object to perform the bilogarithmic methods for determining
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

    def getSigmaP(self, range2fitNCL=None, opt=1):
        """Return the value of the preconsolidation pressure or yield stress.

        Parameters
        ----------
        idx2fitPoly : list, tuple or array (length=2), optional
            Pair of indices between which the 3rd-order polynomial fit is made
            on data to obtain the tangent of the maximum curvature point.
            It works on data regardless of unloading steps.
            The default is [1, 7].
        idx2fitNCL : list, tuple or array (length=2), optional
            Pair of indices between which the 1st-order polynomial fit is made
            on data of normally consolidated line to intersect the bisector
            line and obtain the preconsolidation pressure.
            The default is [7, 11].

        Returns
        -------
        preconsPressure : TYPE
            Preconsolidation pressure or yield stress.
        fig : matplotlib.figure.Figure
            Void ratio at the maximum curvature point.

        """
        sigmaVidx = self.data.findStressIdx(
            stress2find=self.data.sigmaV, cleanedData=True)
        if range2fitNCL is None:  # Indices for fitting the NCL line
            idxInitNCL = -3
            idxEndNCL = None
        else:
            idxInitNCL = self.data.findStressIdx(
                stress2find=range2fitNCL[0], cleanedData=True)
            idxEndNCL = self.data.findStressIdx(
                stress2find=range2fitNCL[1], cleanedData=True) if \
                range2fitNCL[1] < self.data.cleaned['stress'].max() else None

        self.data.raw['specVol'] = self.data.raw['e'] + 1
        self.data.clean()  # -- Updating data without unloads

        def transformX(x, opt=1):
            return np.log(x) if opt == 1 else np.log10(x)

        def transformY(y, opt=1):
            return np.log10(y) if opt == 2 else np.log(y)

        def ticks(x, pos): return f'$e^{np.log(x):.0f}$'

        xScl = 2 if opt == 1 else 10  # log scale for plotting
        yLabel = r'$\log_{10} (1+e)$' if opt == 2 else r'$\ln (1+e)$'

        # -- Linear regresion of points on the preyield line (PYL)
        sigmaPYL = self.data.cleaned['stress'][1: sigmaVidx+1]
        volPYL = self.data.cleaned['specVol'][1: sigmaVidx+1]
        sigmaPYLlog = transformX(sigmaPYL, opt)
        volPYLlog = transformY(volPYL, opt)
        p1_0, p1_1 = polyfit(sigmaPYLlog, volPYLlog, deg=1)
        r2PYL = r2_score(
            y_true=volPYLlog, y_pred=polyval(sigmaPYLlog, [p1_0, p1_1]))
        xPYLFit = np.linspace(
            sigmaPYL.iloc[0], self.data.cleaned['stress'][sigmaVidx+2])
        yPYLFit = polyval(transformX(xPYLFit, opt), [p1_0, p1_1])

        # -- Linear regresion of points on the normally consolidated line (NCL)
        sigmaNCL = self.data.cleaned['stress'][idxInitNCL: idxEndNCL]
        volNCL = self.data.cleaned['specVol'][idxInitNCL: idxEndNCL]
        sigmaNCLLlog = transformX(sigmaNCL, opt)
        volNCLlog = transformY(volNCL, opt)
        p2_0, p2_1 = polyfit(sigmaNCLLlog, volNCLlog, deg=1)
        r2NCL = r2_score(
            y_true=volNCLlog, y_pred=polyval(sigmaNCLLlog, [p2_0, p2_1]))
        xNCLFit = np.linspace(self.data.sigmaV, sigmaNCL.iloc[-1])
        yNCLFit = polyval(transformX(xNCLFit, opt), [p2_0, p2_1])

        # -- Preconsolitadion pressure
        self.sigmaP = xScl ** ((p2_0 - p1_0) / (p1_1 - p2_1))
        vSigmaP = polyval(transformX(self.sigmaP, opt), [p1_0, p1_1])
        vSigmaV = polyval(transformX(self.data.sigmaV, opt), [p1_0, p1_1])

        # -- plot compresibility curve
        fig = plt.figure(figsize=[9, 4.8])
        # ax = fig.add_subplot(111)
        ax = fig.add_axes([0.08, 0.12, 0.55, 0.85])
        # ax.semilogx(self.data.raw['stress'][1:],
        #             transformY(self.data.raw['specVol'][1:], opt),
        #             basex=xScl, ls='--', marker='o', lw=1, c='k', mfc='w',
        #             label='Data')
        ax.plot(self.data.raw['stress'][1:],
                transformY(self.data.raw['specVol'][1:], opt),
                ls='--', marker='o', lw=1, c='k', mfc='w', label='Data')
        ax.set_xscale('log', basex=xScl, subsx=[2, 3, 4, 5, 6, 7, 8, 9])
        if opt == 1:
            ax.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
        methods = ['Butterfield', 'Oikawa', r'Onitsuka \textit{et al.}']
        ax.plot(xPYLFit, yPYLFit, ls='--', c='darkcyan', lw=0.8,
                label=f'Preyield line (R$^2={r2PYL:.3f}$)')
        ax.plot(sigmaPYL, volPYLlog, ls='', marker='.', c='darkcyan')
        ax.plot(xNCLFit, yNCLFit, ls='--', c='crimson', lw=0.8,
                label=f'Postyield line (R$^2={r2NCL:.3f}$)')
        ax.plot(sigmaNCL, volNCLlog, ls='', marker='.', c='crimson')
        ax.plot(self.data.sigmaV, vSigmaV, ls='', marker='|', c='r', ms=15,
                mfc='w', label=str().join([r'$\sigma^\prime_\mathrm{v0}=$ ',
                                           f'{self.data.sigmaV:.0f} kPa']))
        ax.plot(self.sigmaP, vSigmaP, ls='', marker='D', c='r', ms=5,
                mfc='w', label=str().join([r'$\sigma^\prime_\mathrm{p}=$ ',
                                           f'{self.sigmaP:.0f} kPa']))
        # other details
        ax.set(ylabel=f'Specific volume, {yLabel}',
               xlabel=str().join(['Vertical effective stress ',
                                  r'$(\sigma^\prime_\mathrm{v})$ [kPa]']))
        ax.xaxis.set_major_formatter(mtick.ScalarFormatter())
        ax.grid(True, which="both", ls='--', lw=0.5)
        # ax.grid(True, ls='--', lw=0.5)
        ax.legend(bbox_to_anchor=(1.04, 0.5), loc=6, title=str().join([
            r"\textbf{", f"{methods[opt-1]}", "'s method}"]))
        return fig


# %%
'''
BSD 2 license.

Copyright (c) 2020, Exneyder A. Montoya Araque and Alan J. Aparicio-Ortube.
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
'''
