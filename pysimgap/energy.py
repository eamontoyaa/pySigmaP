"""
``energy.py`` module.

Contains the class and its methods for determinig the preconsolidation
pressure from a consolidation test by the energy methods proposed by
Becker et al. (1987), Morin (1988) and Wang & Frost (2004).

References
----------
Becker, D. E., Crooks, J. H. A., Been, K., & Jefferies, M. G. (1987). Work as a
criterion for determining in situ and yield stresses in clays. Canadian
Geotechnical Journal, 24, 4, 549-564, https://doi.org/10.1139/t87-070

Morin, P. (1988). Work as a criterion for determining in situ and yield
stresses in clays: Discussion. Canadian Geotechnical Journal, 25, 4, 845-847,
https://doi.org/10.1139/t88-096

Wang, L. B., & Frost, J. D. (2004). Dissipated strain energy method for
determining preconsolidation pressure. Canadian Geotechnical Journal, 41, 4,
760-768, https://doi.org/10.1139/t04-013

"""

# -- Required modules
import numpy as np
from numpy.polynomial.polynomial import polyfit, polyval
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from sklearn.metrics import r2_score

plt.rcParams['font.family'] = 'Serif'
plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True


class BeckerEtAl():
    """
    ``BeckerEtAl`` class.

    When the object is instanced, the method ``getSigmaP()`` calculates the
    preconsolidation pressure by the method proposed by Becker et al. (1987)
    based on the parameters of the method. See the method documentation for
    more information.

    Attributes
    ----------
    data : Object instanced from the ``Data`` class.
        Contains the data structure from the consolidation test. See the class
        documentation for more information.

    Examples
    --------
    >>> data = Data(pd.read_csv('testData/testData.csv'), sigmaV=75)
    >>> method = BeckerEtAl(data)
    >>> method.getSigmaP(range2fitSCR=None, range2fitLCR=None, zoom=4)
    >>> method.sigmaP, method.ocr
    (670.2104847956236, 8.936139797274983)
    >>> method.getSigmaP(range2fitSCR=None, range2fitLCR=[700, 10000], zoom=4)
    >>> method.sigmaP, method.ocr
    (500.09903645877176, 6.667987152783623)
    """

    def __init__(self, data):
        """Initialize the class."""
        self.data = data
        self.calculateWork()
        return

    def calculateWork(self, morinFormulation=False):
        """
        Calculate the cumulated (total) strain energy.

        Returns
        -------
        None.

        """
        sigma = self.data.raw['stress'].array
        epsilon = self.data.raw['strain'].array
        deltaWork = 0.5*(epsilon[1:] - epsilon[:-1])*(sigma[1:] + sigma[:-1])
        deltaWork = np.hstack((0, deltaWork))
        if morinFormulation:  # Work per unit volume of solids
            deltaWork *= (1 + self.data.e_0)
        self.data.raw['deltaWork'] = deltaWork
        self.data.raw['work'] = np.cumsum(deltaWork)
        self.data.clean()  # Data without unloads
        return

    def getSigmaP(self, range2fitSCR=None, range2fitLCR=None, zoom=3,
                  morinFormulation=False):
        """
        Return the value of the preconsolidation pressure or yield stress.

        Parameters
        ----------
        range2fitSCR : list, tuple or array (length=2), optional
            Initial and final pressures between which the first-order
            polynomial will be fit to the data on the small compressibility
            range (SCR) (i.e. pre-yield range). If None, the SCR will be fit
            from the first point of the curve to the point before the
            in-situ vertical effective stress. The default is None.
        range2fitLCR : list, tuple or array (length=2), optional
            Initial and final pressures between which the first-order
            polynomial will be fit to the data on the large compressibility
            range (LCR) (i.e. post-yield range). If None, the LCR will be
            automatically fit to the same points used for calculating the
            compression index only if it was calculated with a linear fit,
            otherwise, the last three points of the data will be used.
            The default is None.
        zoom : int, optional
            Value to magnify the view of the firsts points of the test and the
            preconsolidation pressure in an inset window. The default is 3.
        morinFormulation : bool, optional
            Boolean to specify if the work per unit volume of solids (Morin
            formulation) must be used instead of the work per unit volume
            (Becker formulation). The default is False.

        Returns
        -------
        fig : matplotlib figure
            Figure with the development of the method and the results.

        """
        # Calculating work
        if morinFormulation:
            self.morinFormulation = morinFormulation
            self.calculateWork()  # Calculate again with Morin Formulation

        # -- Preyield range or small compressibility range (SCR)
        maskSCR = np.full(len(self.data.cleaned), False)
        if range2fitSCR is None:  # Indices for fitting the SCR line
            idxInitSCR = 0
            idxEndSCR = self.data.findStressIdx(
                stress2find=self.data.sigmaV, cleanedData=True)
        else:
            idxInitSCR = self.data.findStressIdx(
                stress2find=range2fitSCR[0], cleanedData=True)
            idxEndSCR = self.data.findStressIdx(
                stress2find=range2fitSCR[1], cleanedData=True)
        maskSCR[idxInitSCR: idxEndSCR] = True
        # -- Linear regresion of points on the preyield line (SCR)
        sigmaSCR = self.data.cleaned['stress'][maskSCR]
        workSCR = self.data.cleaned['work'][maskSCR]
        p1_0, p1_1 = polyfit(sigmaSCR, workSCR, deg=1)
        r2SCR = r2_score(
            y_true=workSCR, y_pred=polyval(sigmaSCR, [p1_0, p1_1]))
        xSCR = np.linspace(0, self.data.cleaned['stress'].iloc[-3])
        ySCR = polyval(xSCR, [p1_0, p1_1])

        # -- Post yield range or large compressibility range
        maskLCR = np.full(len(self.data.cleaned), False)
        if range2fitLCR is not None or self.data.fitCc:
            if range2fitLCR is not None:
                idxInitLCR = self.data.findStressIdx(
                    stress2find=range2fitLCR[0], cleanedData=True)
                idxEndLCR = self.data.findStressIdx(
                    stress2find=range2fitLCR[1], cleanedData=True)
                maskLCR[idxInitLCR: idxEndLCR] = True
            elif self.data.fitCc:
                maskLCR = self.data.maskCc
            # -- Linear regresion of points on post yield line
            sigmaLCR = self.data.cleaned['stress'][maskLCR]
            workLCR = self.data.cleaned['work'][maskLCR]
            lcrInt, lcrSlope = polyfit(sigmaLCR, workLCR, deg=1)
            r2LCR = r2_score(
                y_true=workLCR, y_pred=polyval(sigmaLCR, [lcrInt, lcrSlope]))
        else:  # Using the steepest point of a cubic spline
            sigma = self.data.cleaned['stress']
            cs = CubicSpline(x=sigma, y=self.data.cleaned['work'])
            sigmaCS = np.linspace(sigma.iloc[0], sigma.iloc[-1], 500)
            steepestSlopeIdx = np.argmax(cs(sigmaCS, 1))
            lcrSlope = cs(sigmaCS, 1)[steepestSlopeIdx]
            lcrInt = cs(sigmaCS[steepestSlopeIdx]) - \
                lcrSlope*sigmaCS[steepestSlopeIdx]
        xLCR = np.linspace(
            -lcrInt/lcrSlope, self.data.cleaned['stress'].iloc[-1])
        yLCR = polyval(xLCR, [lcrInt, lcrSlope])

        # -- Preconsolitadion pressure
        workSigmaV = polyval(self.data.sigmaV, [p1_0, p1_1])
        self.sigmaP = (lcrInt - p1_0) / (p1_1 - lcrSlope)
        self.wSigmaP = polyval(self.sigmaP, [p1_0, p1_1])
        self.ocr = self.sigmaP / self.data.sigmaV

        # -- plot compresibility curve
        fig = plt.figure(figsize=[9, 4.8])
        ax = fig.add_axes([0.08, 0.12, 0.65, 0.85])
        ax.plot(self.data.raw['stress'], self.data.raw['work'], ls='--',
                marker='o', lw=1, c='k', mfc='w', label='Data')  # all data

        # Small compressibility range
        ax.plot(xSCR, ySCR, ls='--', c='darkred', lw=0.8,
                label='Preyield line')
        ax.plot(sigmaSCR, workSCR, ls='', marker='+', c='darkred',
                label=f'Data for linear fit\n(R$^2={r2SCR:.3f}$)')
        # Large compressibility range
        ax.plot(xLCR, yLCR, ls='--', c='darkgreen', lw=0.8,
                label='Postyield line')
        if range2fitLCR is not None or self.data.fitCc:
            ax.plot(sigmaLCR, workLCR, ls='', marker='x', c='darkgreen',
                    label=f'Data for linear fit\n(R$^2={r2LCR:.3f}$)')
        # Other plots
        ax.plot(self.data.sigmaV, workSigmaV, ls='', marker='|', c='r', ms=15,
                mfc='w', label=str().join([r'$\sigma^\prime_\mathrm{v0}=$ ',
                                           f'{self.data.sigmaV:.0f} kPa']))
        ax.plot(self.sigmaP, self.wSigmaP, ls='', marker='D', c='r', ms=5,
                mfc='w', label=str().join([r'$\sigma^\prime_\mathrm{p}=$ ',
                                           f'{self.sigmaP:.0f} kPa\n',
                                           f'OCR= {self.ocr:.1f}']))
        # Other details
        methodTitle = r"\textbf{Becker \textit{et al.} method}"
        if morinFormulation:
            methodTitle = r"\textbf{Morin method}"
            yAxisLabel = 'Total work per unit vol. of soilds (W) [kJ m$^{-3}$]'
        else:
            yAxisLabel = 'Total work per unit vol. (W) [kJ m$^{-3}$]'
        ax.set(ylabel=yAxisLabel, xlabel=str().join([
            'Vertical effective stress ',
            r'$(\sigma^\prime_\mathrm{v})$ [kPa]']))
        ax.grid(True, ls='--', lw=0.5)
        ax.legend(bbox_to_anchor=(1.02, 0.5), loc=6, title=methodTitle)

        # -- inset axes to zoom
        axins = zoomed_inset_axes(ax, zoom=zoom, loc=4)
        axins.plot(self.data.raw['stress'], self.data.raw['work'], ls='--',
                   marker='o', lw=1, c='k', mfc='w')  # all data
        axins.plot(xSCR, ySCR, ls='--', c='darkred', lw=0.8)
        axins.plot(sigmaSCR, workSCR, ls='', marker='+', c='darkred')
        axins.plot(xLCR, yLCR, ls='--', c='darkgreen', lw=0.8)
        if range2fitLCR is not None or self.data.fitCc:
            axins.plot(sigmaLCR, workLCR, ls='', marker='x', c='darkgreen')
        axins.plot(self.data.sigmaV, workSigmaV, ls='', marker='|', c='r',
                   ms=15, mfc='w')
        axins.plot(self.sigmaP, self.wSigmaP, marker='D', c='r', ms=5, mfc='w')
        axins.grid(True, ls='--', lw=0.5)
        axins.set(xlim=(-0.05*self.sigmaP, 1.25 * self.sigmaP),
                  ylim=(-0.05 * self.wSigmaP, 2.05 * self.wSigmaP),
                  xlabel=r'$\sigma^\prime_\mathrm{v}$ [kPa]',
                  ylabel=r'W [kJ m$^{-3}$]')
        axins.xaxis.tick_top()
        axins.xaxis.set_label_position('top')
        axins.xaxis.set_tick_params(labelsize='small')
        axins.yaxis.set_tick_params(labelsize='small')
        # axins.set_yticks([])
        mark_inset(ax, axins, loc1=2, loc2=3, fc='none', ec="0.5")  # bbox
        return fig


class WangAndFrost(BeckerEtAl):
    """
    ``WangAndFrost`` class.

    When the object is instanced, the method ``getSigmaP()`` calculates the
    preconsolidation pressure by the method proposed by Wang & Frost (2004)
    based on the parameters of the method. See the method documentation for
    more information.

    Attributes
    ----------
    data : Object instanced from the ``Data`` class.
        Contains the data structure from the consolidation test. See the class
        documentation for more information.

    Examples
    --------
    >>> data = Data(pd.read_csv('testData/testData.csv'), sigmaV=75)
    >>> method = WangAndFrost(data)
    >>> method.getSigmaP(range2fitLCR=None)
    >>> method.sigmaP, method.ocr
    (930.59421331855, 12.407922844247334)
    >>> method.getSigmaP(range2fitLCR=[700, 10000])
    >>> method.sigmaP, method.ocr
    (718.9365936678746, 9.585821248904995)
    """

    def __init__(self, data):
        """Initialize the class."""
        # Heritage BeckerEtAl class
        BeckerEtAl.__init__(self, data)
        return

    def calculateDissipatedE(self, range2fitLCR=None):
        """
        Calculate the accumulative dissipated strain energy corrected.

        Obtain the correction value based on the linear regression on the ADSE
        line or the tangent line to the steepest point of the cubic spline that
        passes through the points.

        Parameters
        ----------
        range2fitLCR : list, tuple or array (length=2), optional
            Initial and final pressures between which the first-order
            polynomial will be fit to the data on the large compressibility
            range (LCR) (i.e. post-yield range). If None, the LCR will be
            automatically obtained with the same criteria used for calculating
            the recompression index with the ``Data`` class. The default is
            None.

        Returns
        -------
        None.

        """
        # -- Incremental total strain energy (ITSE = deltaWork)
        self.data.raw['ITSE'] = self.data.raw['deltaWork']
        # -- Accumulative total strain energy (ATSE = work)
        self.data.raw['ATSE'] = self.data.raw['work']
        # -- Calculating the accumulative elastic strain energy (AESE)
        self.data.raw['AESE'] = (
            self.data.raw['stress'] * self.data.idxCr / (1 + self.data.e_0))
        # -- Calculating the accumulative dissipated strain energy (ADSE)
        self.data.raw['ADSE'] = self.data.raw['ATSE'] - self.data.raw['AESE']
        # -- Calculating value for correcting ADSE (OR: corrVal)
        self.data.clean()  # Updating data without unloads
        if range2fitLCR is not None or self.data.fitCc:
            s2fitTE = self.data.cleaned['stress'][self.maskLCR]
            w2fitTE = self.data.cleaned['ATSE'][self.maskLCR]
            self.corrVal, _ = polyfit(s2fitTE, w2fitTE, deg=1)
        else:
            cs = CubicSpline(x=self.data.cleaned['stress'],
                             y=self.data.cleaned['ATSE'])
            sigmaCS = np.linspace(0, self.data.cleaned['stress'].iloc[-1], 500)
            steepestSlopeIdx = np.argmax(cs(sigmaCS, 1))
            lcrSlope = cs(sigmaCS, 1)[steepestSlopeIdx]
            self.corrVal = cs(sigmaCS[steepestSlopeIdx]) - \
                lcrSlope*sigmaCS[steepestSlopeIdx]
        # -- Calculating the accumulative total strain energy corrected (ATSEC)
        self.data.raw['ATSEC'] = self.data.raw['ATSE'] + abs(self.corrVal)
        # -- Accumulative dissipated strain energy corrected (ADSEC)
        self.data.raw['ADSEC'] = self.data.raw['ADSE'] + abs(self.corrVal)
        # -- Updating data without unloads
        self.data.clean()
        return

    def getSigmaP(self, range2fitLCR=None):
        """
        Return the value of the preconsolidation pressure or yield stress.

        Parameters
        ----------
        range2fitLCR : list, tuple or array (length=2), optional
            Initial and final pressures between which the first-order
            polynomial will be fit to the data on the large compressibility
            range (LCR) (i.e. post-yield range). If None, the LCR will be
            automatically obtained with the same criteria used for calculating
            the recompression index with the ``Data`` class. The default is
            None.

        Returns
        -------
        fig : matplotlib figure
            Figure with the development of the method and the results.

        """
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
            # -- Calculatting the dissipated strain energy corrected (ADSEC)
            self.calculateDissipatedE(range2fitLCR)
            # -- Linear regresion of points on post yield line
            sigmaLCR = self.data.cleaned['stress'][self.maskLCR]
            disspE = self.data.cleaned['ADSEC'][self.maskLCR]
            lcrInt, lcrSlope = polyfit(sigmaLCR, disspE, deg=1)
            r2LCR = r2_score(
                y_true=disspE, y_pred=polyval(sigmaLCR, [lcrInt, lcrSlope]))

        else:  # Using the steepest point of a cubic spline
            # -- Calculatting the dissipated strain energy corrected (ADSEC)
            self.calculateDissipatedE(range2fitLCR)
            sigma = self.data.cleaned['stress']
            cs = CubicSpline(x=sigma, y=self.data.cleaned['ADSEC'])
            sigmaCS = np.linspace(sigma.iloc[0], sigma.iloc[-1], 500)
            steepestSlopeIdx = np.argmax(cs(sigmaCS, 1))
            lcrSlope = cs(sigmaCS, 1)[steepestSlopeIdx]
            lcrInt = cs(sigmaCS[steepestSlopeIdx]) - \
                lcrSlope*sigmaCS[steepestSlopeIdx]

        # -- Preconsolitadion pressure
        self.sigmaP = (abs(self.corrVal) - lcrInt) / lcrSlope
        self.sseSigmaP = polyval(self.sigmaP, [lcrInt, lcrSlope])
        self.ocr = self.sigmaP / self.data.sigmaV

        xLCR = np.linspace(self.sigmaP, self.data.cleaned['stress'].iloc[-1])
        yLCR = polyval(xLCR, [lcrInt, lcrSlope])

        # -- Plot compresibility curve
        fig = plt.figure(figsize=[9, 4.8])
        ax = fig.add_axes([0.08, 0.12, 0.65, 0.85])
        ax.plot(self.data.raw['stress'], self.data.raw['ATSEC'], ls=':',
                marker='v', lw=0.3, c='gray', mfc='w', label='Total energy')
        ax.plot(self.data.raw['stress'], self.data.raw['AESE'], ls=':',
                marker='s', lw=0.3, c='gray', mfc='w', label='Elastic energy')
        ax.plot(self.data.raw['stress'], self.data.raw['ADSEC'], ls='--',
                marker='o', lw=1, c='k', mfc='w', label='Dissipated energy')
        ax.hlines(y=abs(self.corrVal), xmin=0, xmax=self.sigmaP, lw=0.8,
                  color='darkred')
        # Large compressibility range
        ax.plot(xLCR, yLCR, ls='--', c='darkgreen', lw=0.8,
                label='Postyield line')
        if range2fitLCR is not None or self.data.fitCc:
            ax.plot(sigmaLCR, disspE, ls='', marker='x', c='darkgreen',
                    label=f'Data for linear fit\n(R$^2={r2LCR:.3f}$)')
        # Other plots
        ax.plot(self.data.sigmaV, -self.corrVal, ls='', marker='|', c='r',
                ms=15, mfc='w', label=str().join([
                    r'$\sigma^\prime_\mathrm{v0}=$ ',
                    f'{self.data.sigmaV:.0f} kPa']))
        ax.plot(self.sigmaP, self.sseSigmaP, ls='', marker='D', c='r', ms=5,
                mfc='w', label=str().join([r'$\sigma^\prime_\mathrm{p}=$ ',
                                           f'{self.sigmaP:.0f} kPa\n',
                                           f'OCR= {self.ocr:.1f}']))
        # Other details
        ax.set(ylabel='Specific strain energy [kPa]',
               xlabel=str().join(['Vertical effective stress ',
                                  r'$(\sigma^\prime_\mathrm{v})$ [kPa]']))
        ax.grid(True, ls='--', lw=0.5)
        ax.legend(bbox_to_anchor=(1.02, 0.5), loc=6,
                  title=r"\textbf{Wang \& Frost method}")
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
