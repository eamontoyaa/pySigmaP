"""
``casagrande.py`` module.

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
    """``Casagrande`` class.

    When the object is instanced, the method ``getSigmaP()`` calculates the
    preconsolidation pressure by the method proposed by Casagrande (1963)
    based on the compression index obtained with the ``Data`` class and the
    parameters of the method. See the method documentation for more
    information.

    Attributes
    ----------
    data : Object instanced from the ``Data`` class.
        Contains the data structure from the consolidation test. See the class
        documentation for more information.

    Examples
    --------
    >>> data = Data(pd.read_csv('testData/testData.csv'), sigmaV=75)
    >>> method = Casagrande(data)
    >>> method.getSigmaP(range2fitFOP=None, loglog=True)
    >>> method.sigmaP, method.ocr
    (925.6233444375316, 12.34164459250042)
    >>> method.getSigmaP(range2fitFOP=[20, 5000], loglog=True)
    >>> method.sigmaP, method.ocr
    (651.5675456910274, 8.687567275880365)
    >>> method.getSigmaP(mcp=200)
    >>> method.sigmaP, method.ocr
    (498.6050168481866, 6.6480668913091545)
    """

    def __init__(self, data):
        """Initialize the class."""
        self.data = data
        return

    def getSigmaP(self, mcp=None, range2fitFOP=None, loglog=True):
        """
        Return the value of the preconsolidation pressure or yield stress.

        Parameters
        ----------
        mcp : float, optional
            Variable to manually specify the maximun curvature point. If a
            value is input, the other parameters are not taken into account.
            The default is None.
        range2fitFOP : list, tuple or array (length=2), optional
            Initial and final pressures between which the fourth-order
            polynomial (FOP) will be fit to the compressibility curve to
            calculate the maximum curvature point (MCP). If None, the MCP will
            be obtained using a cubic spline that passes through the data. The
            default is None.
        loglog : bool, optional
            Boolean to specify if the vertical effective stress will be
            transformed applying a logarithm twice to fit the FOT. The default
            is True.

        Returns
        -------
        fig : matplotlib figure
            Figure with the development of the method and the results.

        """
        def transform(x, reverse=False):
            if reverse:  # Remove a logaritmic scale
                if loglog:
                    return 10**10**x
                else:
                    return 10**x
            else:  # Set a logaritmic scale
                if loglog:
                    return np.log10(np.log10(x))
                else:
                    return np.log10(x)

        sigmaLog = np.log10(self.data.cleaned['stress'][1:])
        cs = CubicSpline(x=sigmaLog, y=self.data.cleaned['e'][1:])
        sigmaCS = np.linspace(sigmaLog.iloc[0], sigmaLog.iloc[-1], 100)
        if range2fitFOP is None:  # Using a cubic spline
            x4FOP = 10**sigmaCS
            # -- Curvature function k(x) = f''(x)/(1+(f'(X))²)³/²
            curvature = abs(cs(sigmaCS, 2)) / (1+cs(sigmaCS, 1)**2)**(3/2)
            maxCurvIdx = find_peaks(curvature, distance=500)[0][0]
            # maxCurvIdx = np.argmax(curvature)  # Max. Curvature index
            self.sigmaMC = 10**sigmaCS[maxCurvIdx]  # Max. Curvature point
            self.eMC = cs(sigmaCS[maxCurvIdx])  # Void ratio at MC

        else:  # Using a fourth order polynomial (FOP)
            # -- Indices to fit the FOP
            idxInitFOP = self.data.findStressIdx(
                stress2find=range2fitFOP[0], cleanedData=True)
            idxEndFOP = self.data.findStressIdx(
                stress2find=range2fitFOP[1], cleanedData=True)

            # -- fittig a polynomial to data without unloads
            sigmaFOP = self.data.cleaned['stress'][idxInitFOP: idxEndFOP]
            sigmaFOPlog = transform(sigmaFOP)
            eFOP = self.data.cleaned['e'][idxInitFOP: idxEndFOP]
            p0, p1, p2, p3, p4 = polyfit(sigmaFOPlog, eFOP, deg=4)
            r2FOP = r2_score(y_true=eFOP,
                             y_pred=polyval(sigmaFOPlog, [p0, p1, p2, p3, p4]))
            x4FOP = np.linspace(sigmaFOP.iloc[0], sigmaFOP.iloc[-1], 1000)
            x4FOPlog = transform(x4FOP)
            y4FOP = polyval(x4FOPlog, [p0, p1, p2, p3, p4])

            # -- Curvature function k(x) = f''(x)/(1+(f'(X))²)³/²
            firstDer = p1 + 2*p2*x4FOPlog + 3*p3*x4FOPlog**2 + 4*p4*x4FOPlog**3
            secondDer = 2*p2 + 6*p3*x4FOPlog + 12*p4*x4FOPlog**2
            curvature = abs(secondDer) / (1+firstDer**2)**1.5
            maxCurvIdx = np.argmax(curvature)    # Max. Curvature index
            self.sigmaMC = transform(x4FOPlog[maxCurvIdx], True)  # Max. Curv.
            self.eMC = y4FOP[maxCurvIdx]  # Void ratio at max. curvature

        if mcp is not None:
            self.sigmaMC = mcp  # Max. Curvature point
            self.eMC = cs(np.log10(self.sigmaMC))  # Void ratio at MC

        # -- Bisector line
        slopeMP = cs(np.log10(self.sigmaMC), nu=1)  # Slope at MC
        y1, x1 = self.eMC, np.log10(self.sigmaMC)
        x2 = np.log10(np.linspace(
            self.sigmaMC, self.data.cleaned['stress'].iloc[-1], 50))
        y2 = slopeMP * (x2 - x1) + y1
        slopeBisect = np.tan(0.5*np.arctan(slopeMP))  # slope of bisector line
        y2bis = slopeBisect * (x2 - x1) + y1

        # -- Preconsolidation pressure
        self.sigmaP = 10**((y1 - slopeBisect*x1 - self.data.idxCcInt) /
                           (-self.data.idxCc - slopeBisect))
        self.eSigmaP = slopeBisect * (np.log10(self.sigmaP) - x1) + y1
        self.ocr = self.sigmaP / self.data.sigmaV

        # -- plotting
        fig = plt.figure(figsize=[9, 4.8])
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
        # Compression index (Cc)
        x4Cc = np.linspace(self.sigmaP, self.data.cleaned['stress'].iloc[-1])
        y4Cc = -self.data.idxCc * np.log10(x4Cc) + self.data.idxCcInt
        l2 = ax1.plot(x4Cc, y4Cc, ls='--', lw=0.8, color='darkgreen',
                      label=str().join([r'$C_\mathrm{c}=$',
                                        f'{self.data.idxCc:.3f}']))
        # Curvature
        l5 = ax2.plot(x4FOP, curvature, ls='--', c='darkorange', lw=0.8,
                      mfc='w', label='Curvature')
        # Other plots
        l6 = ax1.plot(self.sigmaMC, self.eMC, ls='', marker='o', c='r',
                      mfc='w', label='Max. curvature point')
        l7 = ax1.plot(self.sigmaP, self.eSigmaP, ls='', marker='D', c='r',
                      ms=5, mfc='w',
                      label=str().join([r'$\sigma^\prime_\mathrm{p}=$ ',
                                        f'{self.sigmaP:.0f} kPa\n',
                                        f'OCR= {self.ocr:.1f}']))
        # Fourth order polynomial fit
        if range2fitFOP is not None:
            l3 = ax1.plot(x4FOP, y4FOP, ls='--', lw=0.8, color='darkred',
                          label=r'$4^\mathrm{th}$-order polynomial')
            l4 = ax1.plot(sigmaFOP, eFOP, ls='', marker='x', c='darkred',
                          label=f'Data for linear fit\n(R$^2={r2FOP:.3f}$)')
            allLayers = l1 + l2 + l3 + l4 + l5 + l6 + l7
        else:
            allLayers = l1 + l2 + l5 + l6 + l7
        # Legend
        labs = [layer.get_label() for layer in allLayers]
        ax2.legend(allLayers, labs, bbox_to_anchor=(1.15, 0.5), loc=6,
                   title=r"\textbf{Casagrande method}")
        # Other details
        ax1.set(xscale='log', ylabel=r'Void ratio $(e)$',
                xlabel=str().join(['Vertical effective stress ',
                                  r'$(\sigma^\prime_\mathrm{v})$ [kPa]']))
        ax1.xaxis.set_major_formatter(mtick.ScalarFormatter())
        ax2.set(ylabel='Curvature $(k)$')
        ax1.grid(True, which="both", ls='--', lw=0.5)
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
