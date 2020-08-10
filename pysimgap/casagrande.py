"""
casagrande.py module.

Contains the class and its methods for determinig the preconsolidation
pressure from a consolidation test by the method proposed by Casagrande (1963).

References
----------
Casagrande, A. (1936). The determination of pre-consolidation load and its
practical significance. In Proceedings of the First International Conference
on Soil Mechanins and Foundations Engineering, 3, 60-64.

"""

# -- Required modules
import numpy as np
from numpy.polynomial.polynomial import polyfit, polyval
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import matplotlib.ticker as mtick

plt.rcParams['font.family'] = 'Serif'
plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True


class Casagrande:
    """Casagrande class."""

    def __init__(self, data):
        """
        Initialize the Casagrande class.

        Instance an object to perform the Casagrande's method for determining
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
        Return the value of the preconsolidation pressure or yield stress.

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
        if None in [range2fitTOP, range2fitNCL]:  # Using a cubic spline
            sigmaLog = np.log10(self.data.cleaned['stress'][1:])
            cs = CubicSpline(x=sigmaLog, y=self.data.cleaned['e'][1:])
            sigmaCS = np.linspace(sigmaLog.iloc[0], sigmaLog.iloc[-1], 500)
            xFitTOP = 10**sigmaCS  # Just to plot the spline

            # -- Curvature function k(x) = f''(x)/(1+(f'(X))²)³/²
            curvature = abs(cs(sigmaCS, 2)) / (1+cs(sigmaCS, 1)**2)**(3/2)
            maxCurvIdx = find_peaks(curvature, distance=500)[0][0]
            self.sigmaMC = 10**sigmaCS[maxCurvIdx]  # Max. Curvature ppint (MC)
            self.eMC = cs(sigmaCS[maxCurvIdx])  # Void ratio at MC
            pPrime = cs(np.log10(self.sigmaMC), nu=1)  # Slope at MC

            # NCL line
            steepestSlopeIdx = np.argmin(cs(sigmaCS, 1))
            slopeNCL = cs(sigmaCS, 1)[steepestSlopeIdx]
            intNCL = cs(sigmaCS[steepestSlopeIdx]) - \
                slopeNCL*sigmaCS[steepestSlopeIdx]

        if range2fitTOP is not None:  # Using a third order poly (TOP) fit
            # -- Indices to fit the TOP
            idxInitTOP = self.data.findStressIdx(
                stress2find=range2fitTOP[0], cleanedData=True)
            idxEndTOP = self.data.findStressIdx(
                stress2find=range2fitTOP[1], cleanedData=True)

            # -- fittig a polynomial to data without unloads
            sigmaTOP = self.data.cleaned['stress'][idxInitTOP: idxEndTOP]
            sigmaTOPlog = np.log10(sigmaTOP)
            eTOP = self.data.cleaned['e'][idxInitTOP: idxEndTOP]
            p1_0, p1_1, p1_2, p1_3 = polyfit(sigmaTOPlog, eTOP, deg=3)
            r2TOP = r2_score(y_true=eTOP, y_pred=polyval(
                sigmaTOPlog, [p1_0, p1_1, p1_2, p1_3]))
            xFitTOP = np.linspace(sigmaTOP.iloc[0], sigmaTOP.iloc[-1], 500)
            yFitTOP = polyval(np.log10(xFitTOP), [p1_0, p1_1, p1_2, p1_3])

            # -- Curvature function k(x) = f''(x)/(1+(f'(X))²)³/²
            p2_0, p2_1, p2_2, p2_3 = polyfit(
                np.log10(sigmaTOPlog), eTOP, deg=3)
            xFitTOPloglog = np.log10(np.log10(xFitTOP))
            num = 2*p2_2 + 6*p2_3*xFitTOPloglog
            den = 1 + (p2_1 + 2*p2_2*xFitTOPloglog +
                       3*p2_3*xFitTOPloglog**2)**2
            curvature = abs(num) / den**(3/2)
            maxCurvIdx = list(curvature).index(curvature.max())
            self.sigmaMC = 10**10**xFitTOPloglog[maxCurvIdx]  # Max. Curvature
            self.eMC = yFitTOP[maxCurvIdx]  # Void ratio at max. curvature

            # -- Slope at the maximum curvature point and tangent line
            pPrime = p1_1 + 2*p1_2*np.log10(self.sigmaMC) +\
                3*p1_3*np.log10(self.sigmaMC)**2

        if range2fitNCL is not None:
            # -- Indices to fit the NCL line
            idxInitNCL = self.data.findStressIdx(
                stress2find=range2fitNCL[0], cleanedData=True)
            idxEndNCL = self.data.findStressIdx(
                stress2find=range2fitNCL[1], cleanedData=True) if \
                range2fitNCL[1] < self.data.cleaned['stress'].max() else None

            # -- Linear regresion of points on normally consolidated line (NCL)
            sigmaNCL = self.data.cleaned['stress'][idxInitNCL: idxEndNCL]
            sigmaNCLlog = np.log10(sigmaNCL)
            eNCL = self.data.cleaned['e'][idxInitNCL: idxEndNCL]
            intNCL, slopeNCL = polyfit(sigmaNCLlog, eNCL, deg=1)
            r2NCL = r2_score(
                y_true=eNCL, y_pred=polyval(sigmaNCLlog, [intNCL, slopeNCL]))

        # -- Bisector line
        y1, x1 = self.eMC, np.log10(self.sigmaMC)
        x2 = np.log10(np.linspace(
            self.sigmaMC, self.data.cleaned['stress'].iloc[-1], 500))
        y2 = pPrime * (x2 - x1) + y1
        bisPrime = np.tan(0.5*np.arctan(pPrime))  # slope of bisector line
        y2bis = bisPrime * (x2 - x1) + y1

        # -- NC line
        xFitNCL = np.linspace(
            self.data.sigmaV, self.data.cleaned['stress'].iloc[-1], 500)
        yFitNCL = polyval(np.log10(xFitNCL), [intNCL, slopeNCL])

        # -- Preconsolidation pressure
        self.sigmaP = 10**((y1 - bisPrime*x1 - intNCL) / (slopeNCL - bisPrime))
        self.eSigmaP = bisPrime * (np.log10(self.sigmaP) - x1) + y1

        # -- plotting
        fig = plt.figure(figsize=[9, 4.8])
        # ax1 = fig.add_subplot(111)  # Compressibility curve
        ax1 = fig.add_axes([0.08, 0.12, 0.55, 0.85])
        ax2 = ax1.twinx()  # second y axis for curvature function
        l1 = ax1.plot(self.data.raw['stress'][1:], self.data.raw['e'][1:],
                      ls='--', marker='o', lw=1, c='k', mfc='w',
                      label='Compressibility curve')
        # Lines of the Casagrande's method
        ax1.plot([self.sigmaMC, self.data.cleaned['stress'].iloc[-1]],
                 [self.eMC, self.eMC], ls='-.', lw=0.5, c='k')  # hztl line
        ax1.plot(10**x2, y2, ls='-.', lw=0.5, color='k')  # tangent line
        ax1.plot(10**x2, y2bis, ls='-.', lw=0.5, color='k')  # bisector line
        if range2fitNCL is None:
            l2 = ax1.plot(xFitNCL, yFitNCL, ls='--', lw=0.8, color='crimson',
                          label='Normally consolidated line')
        else:
            ax1.plot(sigmaNCL, eNCL, ls='', marker='.', lw=0.8,
                     color='crimson')
            l2 = ax1.plot(xFitNCL, yFitNCL, ls='--', lw=0.8, color='crimson',
                          label=f'NCL linear fit\n(R$^2={r2NCL:.3f}$)')
        if range2fitTOP is None:
            l3 = ax1.plot(10**sigmaCS, cs(sigmaCS), ls='--', lw=0.8,
                          color='darkcyan', label='Cubic spline')
        else:
            ax1.plot(sigmaTOP, eTOP, ls='', marker='.', lw=0.8,
                     color='darkcyan')
            l3 = ax1.plot(xFitTOP, yFitTOP, ls='--', lw=0.8, color='darkcyan',
                          label=str().join([r'$3^\mathrm{rd}$-order poly',
                                            f' fit \n(R$^2={r2TOP:.3f}$)']))
        l4 = ax2.plot(xFitTOP, curvature, ls='--', c='darkorange', lw=0.8,
                      mfc='w', label='Curvature')
        l5 = ax1.plot(self.sigmaMC, self.eMC, ls='', marker='o', c='r',
                      mfc='w', label='Max. curvature point')
        l6 = ax1.plot(self.sigmaP, self.eSigmaP, ls='', marker='D', c='r',
                      ms=5, mfc='w',
                      label=''.join([r'$\sigma^\prime_\mathrm{p}=$ ',
                                     f'{self.sigmaP:.0f} kPa']))
        # other details
        ax1.set(xscale='log', ylabel=r'Void ratio $(e)$',
                xlabel=str().join(['Vertical effective stress ',
                                  r'$(\sigma^\prime_\mathrm{v})$ [kPa]']))
        ax1.xaxis.set_major_formatter(mtick.ScalarFormatter())
        ax2.set(ylabel='Curvature $(k)$')
        ax1.grid(True, which="both", ls='--', lw=0.5)
        # Legend
        allLayers = l1 + l2 + l3 + l4 + l5 + l6
        labs = [layer.get_label() for layer in allLayers]
        ax2.legend(allLayers, labs, bbox_to_anchor=(1.15, 0.5), loc=6,
                   title=r"\textbf{Casagrande's method}")
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
