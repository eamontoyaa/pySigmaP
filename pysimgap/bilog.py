"""
``bilog.py`` module.

Contains the class and its methods for determinig the preconsolidation
pressure from a consolidation test by the bilogarithminc methods proposed by
Butterfield (1979), Oikawa (1987) and Onitsuka et al. (1995).

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
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import matplotlib.ticker as mtick

plt.rcParams['font.family'] = 'Serif'
plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True


class Bilog():
    """
    ``Bilog`` class.

    When the object is instanced, the method ``getSigmaP()`` calculates the
    preconsolidation pressure by the method proposed by Boone (2010) based on
    the parameters of the method.  See the method documentation for more
    information.

    Attributes
    ----------
    data : Object instanced from the ``Data`` class
        Contains the data structure from the consolidation test. See the class
        documentation for more information.

    Examples
    --------
    >>> data = Data(pd.read_csv('testData/testData.csv'), sigmaV=75)
    >>> method = Bilog(data)
    >>> method.getSigmaP(opt=1)  # Butterfield (1979)
    >>> method.sigmaP, method.ocr
    (433.29446692547555, 5.777259559006341)
    >>> method.getSigmaP(range2fitLCR=[1000, 9000], opt=1) # Butterfield
    >>> method.sigmaP, method.ocr
    (399.3361458736937, 5.324481944982582)
    >>> method.getSigmaP(opt=2)  # Oikawa (1987)
    >>> method.sigmaP, method.ocr
    (433.2944669254731, 5.777259559006308)
    >>> method.getSigmaP(range2fitLCR=[1000, 9000], opt=2)  # Oikawa
    >>> method.sigmaP, method.ocr
    (399.3361458736935, 5.32448194498258)
    >>> method.getSigmaP(opt=3)  # Onitsuka et al. (1995)
    >>> method.sigmaP, method.ocr
    (433.2944669254753, 5.777259559006338)
    >>> method.getSigmaP(range2fitLCR=[1000, 9000], opt=3)  # Onitsuka et al.
    >>> method.sigmaP, method.ocr
    (399.3361458736939, 5.324481944982585)
    """

    def __init__(self, data):
        """Initialize the class."""
        self.data = data
        return

    def getSigmaP(self, range2fitSCR=None, range2fitLCR=None, opt=1):
        """
        Return the value of the preconsolidation pressure or yield stress.

        Parameters
        ----------
        range2fitSCR : list, tuple or array (length=2), optional
            Initial and final pressures between which the first-order
            polynomial will be fit to the compressibility curve on the small
            compressibility range (SCR). If None, the SCR will be fit from the
            first point of the compresibility curve to the point before the
            in-situ vertical effective stress. The default is None.
        range2fitLCR : list, tuple or array (length=2), optional
            Initial and final pressures between which the first-order
            polynomial will be fit to the compressibility curve on the normally
            consolidated line (NCL), also knonw as the "large compressibility
            range (LCR)". If None, the NCL will be automatically fit to the
            last three points of the data.
        opt : int, optional
            Integer value to indicate which bilogarithmic method will be used.
            Butterfield (opt=1), Oikawa (opt=2) and Onitsuka et al. (opt=3).
            The default is 1.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Figure with the development of the method and the results.

        """
        def transformX(x, opt=1):
            return np.log(x) if opt == 1 else np.log10(x)

        def transformY(y, opt=1):
            return np.log10(y) if opt == 2 else np.log(y)

        self.data.raw['vol'] = self.data.raw['e'] + 1
        self.data.clean()  # -- Updating data without unloads

        # def ticks(x, pos): return f'$e^{np.log(x):.0f}$'

        if range2fitSCR is None:  # Indices for fitting the SCR line
            idxInitSCR = 1
            idxEndSCR = self.data.findStressIdx(
                stress2find=self.data.sigmaV, cleanedData=True) - 1
        else:
            idxInitSCR = self.data.findStressIdx(
                stress2find=range2fitSCR[0], cleanedData=True)
            idxEndSCR = self.data.findStressIdx(
                stress2find=range2fitSCR[1], cleanedData=True)

        # -- Linear regresion of points on the small compressibility range(SCR)
        sigmaSCR = self.data.cleaned['stress'][idxInitSCR: idxEndSCR+1]
        volSCR = self.data.cleaned['vol'][idxInitSCR: idxEndSCR+1]
        sigmaSCRlog = transformX(sigmaSCR, opt)
        volSCRlog = transformY(volSCR, opt)
        p1_0, p1_1 = polyfit(sigmaSCRlog, volSCRlog, deg=1)
        r2SCR = r2_score(
            y_true=volSCRlog, y_pred=polyval(sigmaSCRlog, [p1_0, p1_1]))

        if range2fitLCR is None:  # Indices for fitting the NCL line
            idxInitNCL = -3
            idxEndNCL = None
        else:
            idxInitNCL = self.data.findStressIdx(
                stress2find=range2fitLCR[0], cleanedData=True)
            idxEndNCL = self.data.findStressIdx(
                stress2find=range2fitLCR[1], cleanedData=True)

        # -- Cubic spline that passes through the data
        sigmaLog = transformX(self.data.cleaned['stress'][1:], opt)
        volLog = transformY(self.data.cleaned['vol'][1:], opt)
        cs = CubicSpline(x=sigmaLog, y=volLog)
        # Specific volume at sigma V
        volSigmaV = cs(transformX(self.data.sigmaV, opt))

        # -- Post yield range or large compressibility range
        self.maskLCR = np.full(len(self.data.cleaned), False)
        if range2fitLCR is not None or self.data.fitCc:
            if range2fitLCR is not None:
                idxInitLCR = self.data.findStressIdx(
                    stress2find=range2fitLCR[0], cleanedData=True)
                idxEndLCR = self.data.findStressIdx(
                    stress2find=range2fitLCR[1], cleanedData=True)
                self.maskLCR[idxInitLCR: idxEndLCR] = True
            elif self.data.fitCc:
                self.maskLCR = self.data.maskCc
            # -- Linear regresion of points on post yield line
            sigmaLCR = self.data.cleaned['stress'][self.maskLCR]
            sigmaLCRlog = transformX(sigmaLCR, opt)
            volLCR = self.data.cleaned['vol'][self.maskLCR]
            volLCRlog = transformY(volLCR, opt)
            lcrInt, lcrSlope = polyfit(sigmaLCRlog, volLCRlog, deg=1)
            r2LCR = r2_score(y_true=volLCRlog,
                             y_pred=polyval(volLCRlog, [lcrInt, lcrSlope]))

        else:  # Using the steepest point of a cubic spline
            sigmaCS = np.linspace(sigmaLog.iloc[0], sigmaLog.iloc[-1], 500)
            steepestSlopeIdx = np.argmin(cs(sigmaCS, 1))
            lcrSlope = cs(sigmaCS, 1)[steepestSlopeIdx]
            lcrInt = cs(sigmaCS[steepestSlopeIdx]) - \
                lcrSlope*sigmaCS[steepestSlopeIdx]

        xScl = np.e if opt == 1 else 10  # log scale for plotting
        yLabel = r'$\log_{10} (1+e)$' if opt == 2 else r'$\ln (1+e)$'

        # -- Preconsolitadion pressure
        self.sigmaP = xScl ** ((lcrInt - p1_0) / (p1_1 - lcrSlope))
        self.vSigmaP = polyval(transformX(self.sigmaP, opt), [p1_0, p1_1])
        self.ocr = self.sigmaP / self.data.sigmaV

        # -- Lines of the bilogarithmic methods
        xLCR = np.linspace(self.sigmaP, self.data.cleaned['stress'].iloc[-1])
        yLCR = polyval(transformX(xLCR, opt), [lcrInt, lcrSlope])

        # Small compressibility range
        xSCRFit = np.linspace(sigmaSCR.iloc[0], self.sigmaP)
        ySCRFit = polyval(transformX(xSCRFit, opt), [p1_0, p1_1])

        # -- plot compresibility curve
        fig = plt.figure(figsize=[9, 4.8])
        ax = fig.add_axes([0.08, 0.12, 0.55, 0.85])
        ax.semilogx(self.data.raw['stress'][1:],
                    transformY(self.data.raw['vol'][1:], opt),
                    basex=xScl, ls='--', marker='o', lw=1, c='k', mfc='w',
                    label='Data')
        methods = ['Butterfield', 'Oikawa', r'Onitsuka \textit{et al.}']
        # Small compressibility range
        ax.plot(xSCRFit, ySCRFit, ls='--', c='darkred', lw=0.8,
                label='Small compressibility range')
        ax.plot(sigmaSCR, volSCRlog, ls='', marker='+', c='darkred',
                label=f'Data for linear fit\n(R$^2={r2SCR:.3f}$)')
        # Large compressibility range
        ax.plot(xLCR, yLCR, ls='--', c='darkgreen', lw=0.8,
                label='Large compressibility range')
        if range2fitLCR is not None or self.data.fitCc:
            ax.plot(sigmaLCR, volLCRlog, ls='', marker='x', c='darkgreen',
                    label=f'Data for linear fit\n(R$^2={r2LCR:.3f}$)')
        # Other plots
        ax.plot(self.data.sigmaV, volSigmaV, ls='', marker='|', c='r', ms=15,
                mfc='w', label=str().join([r'$\sigma^\prime_\mathrm{v0}=$ ',
                                           f'{self.data.sigmaV:.0f} kPa']))
        ax.plot(self.sigmaP, self.vSigmaP, ls='', marker='D', c='r', ms=5,
                mfc='w', label=str().join([r'$\sigma^\prime_\mathrm{p}=$ ',
                                           f'{self.sigmaP:.0f} kPa\n',
                                           f'OCR= {self.ocr:.1f}']))
        # Other details
        ax.set(ylabel=f'Specific volume, {yLabel}',
               xlabel=str().join(['Vertical effective stress ',
                                  r'$(\sigma^\prime_\mathrm{v})$ [kPa]']))
        ax.xaxis.set_major_formatter(mtick.ScalarFormatter())
        ax.grid(True, which="both", ls='--', lw=0.5)
        ax.legend(bbox_to_anchor=(1.04, 0.5), loc=6, title=str().join([
            r"\textbf{", f"{methods[opt-1]}", "'s method}"]))
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
