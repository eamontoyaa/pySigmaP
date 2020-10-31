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
import matplotlib.ticker as mtick

plt.rcParams['font.family'] = 'Serif'
plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True
# High-contrast qualitative colour scheme
colors = ('#DDAA33',  # yellow
          '#BB5566',  # red
          '#004488')  # blue


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
    >>> urlCSV = ''.join(['https://raw.githubusercontent.com/eamontoyaa/',
    >>>                   'data4testing/main/pysigmap/testData.csv'])
    >>> data = Data(pd.read_csv(urlCSV), sigmaV=75)
    >>> method = BeckerEtAl(data)
    >>> method.getSigmaP(range2fitRR=None, range2fitCR=None, zoom=4)
    >>> method.sigmaP, method.ocr
    (670.2104847956236, 8.936139797274983)
    >>> method.getSigmaP(range2fitRR=None, range2fitCR=[700, 10000], zoom=4)
    >>> method.sigmaP, method.ocr
    (500.09903645877176, 6.667987152783623)
    """

    def __init__(self, data):
        """Initialize the class."""
        self.data = data
        self.morinFormulation = False
        self.calculateWork()
        return

    def calculateWork(self):
        """
        Calculate the work per unit volume.

        Returns
        -------
        None.

        """
        sigma = self.data.raw['stress'].array
        epsilon = self.data.raw['strain'].array
        deltaWork = 0.5*(epsilon[1:] - epsilon[:-1])*(sigma[1:] + sigma[:-1])
        deltaWork = np.hstack((0, deltaWork))
        if self.morinFormulation:  # Work per unit volume of solids
            deltaWork /= (1 + self.data.e_0)
        self.data.raw['deltaWork'] = deltaWork
        self.data.raw['work'] = np.cumsum(deltaWork)
        self.data.clean()  # Data without unloads
        return

    def getSigmaP(self, range2fitRR=None, range2fitCR=None, zoom=3,
                  morinFormulation=False):
        """
        Return the value of the preconsolidation pressure or yield stress.

        Parameters
        ----------
        range2fitRR : list, tuple or array (length=2), optional
            Initial and final pressures between which the first order
            polynomial will be fit to the data on the recompression range (RR)
            (range before the preconsolidation pressure). If None, the first
            order polynomial will be fit from the first point of the curve to
            the point before the in-situ vertical effective stress. The default
            is None.
        range2fitCR : list, tuple or array (length=2), optional
            Initial and final pressures between which the first order
            polynomial will be fit to the data on the compression range (CR)
            (range beyond the preconsolidation pressure). If None, the CR will
            be automatically fit to the same points used for calculating the
            compression index with the ``Data`` class only if it was calculated
            with a linear fit, otherwise, the steepest slope of the cubic
            spline that passes through the data will be used. The default is
            None.
        zoom : int, optional
            Value to magnify the view of the firsts points of the curve and the
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
        # -- Calculating work
        if morinFormulation:
            self.morinFormulation = morinFormulation
            self.calculateWork()  # Calculate again with Morin Formulation

        # -- Preyield range or recompression range (RR)
        maskRR = np.full(len(self.data.cleaned), False)
        if range2fitRR is None:  # Indices for fitting the RR line
            idxInitRR = 0
            idxEndRR = self.data.findStressIdx(
                stress2find=self.data.sigmaV, cleanedData=True)
        else:
            idxInitRR = self.data.findStressIdx(
                stress2find=range2fitRR[0], cleanedData=True)
            idxEndRR = self.data.findStressIdx(
                stress2find=range2fitRR[1], cleanedData=True)
        maskRR[idxInitRR: idxEndRR] = True
        # Linear regresion
        sigmaRR = self.data.cleaned['stress'][maskRR]
        workRR = self.data.cleaned['work'][maskRR]
        p1_0, p1_1 = polyfit(sigmaRR, workRR, deg=1)
        r2RR = r2_score(
            y_true=workRR, y_pred=polyval(sigmaRR, [p1_0, p1_1]))
        xRR = np.linspace(0, self.data.cleaned['stress'].iloc[-3])
        yRR = polyval(xRR, [p1_0, p1_1])

        # -- Post yield range or compression range (CR)
        maskCR = np.full(len(self.data.cleaned), False)
        if range2fitCR is not None or self.data.fitCc:  # Using a linear fit
            if range2fitCR is not None:
                idxInitCR = self.data.findStressIdx(
                    stress2find=range2fitCR[0], cleanedData=True)
                idxEndCR = self.data.findStressIdx(
                    stress2find=range2fitCR[1], cleanedData=True)
                maskCR[idxInitCR: idxEndCR] = True
            elif self.data.fitCc:
                maskCR = self.data.maskCc
            # Linear regresion
            sigmaCR = self.data.cleaned['stress'][maskCR]
            workCR = self.data.cleaned['work'][maskCR]
            lcrInt, lcrSlope = polyfit(sigmaCR, workCR, deg=1)
            r2CR = r2_score(
                y_true=workCR, y_pred=polyval(sigmaCR, [lcrInt, lcrSlope]))
        else:  # Using the steepest point of a cubic spline
            sigma = self.data.cleaned['stress']
            cs = CubicSpline(x=sigma, y=self.data.cleaned['work'])
            sigmaCS = np.linspace(sigma.iloc[0], sigma.iloc[-1], 500)
            steepestSlopeIdx = np.argmax(cs(sigmaCS, 1))
            lcrSlope = cs(sigmaCS, 1)[steepestSlopeIdx]
            lcrInt = cs(sigmaCS[steepestSlopeIdx]) - \
                lcrSlope*sigmaCS[steepestSlopeIdx]
        xCR = np.linspace(
            -lcrInt/lcrSlope, self.data.cleaned['stress'].iloc[-1])
        yCR = polyval(xCR, [lcrInt, lcrSlope])

        # -- Preconsolitadion pressure
        workSigmaV = polyval(self.data.sigmaV, [p1_0, p1_1])
        self.sigmaP = (lcrInt - p1_0) / (p1_1 - lcrSlope)
        self.wSigmaP = polyval(self.sigmaP, [p1_0, p1_1])
        self.ocr = self.sigmaP / self.data.sigmaV

        # -- Plot compresibility curve
        fig = plt.figure(figsize=[9, 4.8])
        ax = fig.add_axes([0.08, 0.12, 0.55, 0.85])
        ax.plot(self.data.raw['stress'], self.data.raw['work'], ls=(0, (1, 1)),
                marker='o', lw=1.5, c='k', mfc='w', label='Data')  # all data
        # Recompression range
        ax.plot(xRR, yRR, ls='-', c=colors[2], lw=1.125,
                label='Recompression range')
        ax.plot(sigmaRR, workRR, ls='', marker='+', c=colors[2],
                label=f'Data for linear fit\n(R$^2={r2RR:.3f}$)')
        # Compression range
        ax.plot(xCR, yCR, ls='-', c=colors[1], lw=1.125,
                label='Compression range')
        if range2fitCR is not None or self.data.fitCc:
            ax.plot(sigmaCR, workCR, ls='', marker='x', c=colors[1],
                    label=f'Data for linear fit\n(R$^2={r2CR:.3f}$)')
        # Other plots
        ax.plot(self.data.sigmaV, workSigmaV, ls='', marker='|', c='r', ms=15,
                mfc='w', mew=1.5,
                label=str().join([r'$\sigma^\prime_\mathrm{v0}=$ ',
                                  f'{self.data.sigmaV:.0f} kPa']))
        ax.plot(self.sigmaP, self.wSigmaP, ls='', c=colors[0], ms=7, mfc='w',
                marker='o', mew=1.5,
                label=str().join([r'$\sigma^\prime_\mathrm{p}=$ ',
                                  f'{self.sigmaP:.0f} kPa\n',
                                  f'OCR= {self.ocr:.1f}']))
        # Other details
        if morinFormulation:
            methodTitle = r"\textbf{Morin method}"
        else:
            methodTitle = r"\textbf{Becker \textit{et al.} method}"
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set(ylabel='Total work per unit volume, $W$ [kJ m$^{-3}$]',
               xlabel=str().join(['Vertical effective stress, ',
                                  r'$\sigma^\prime_\mathrm{v}$ [kPa]']))
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.grid(False)
        ax.legend(bbox_to_anchor=(1.125, 0.5), loc=6, title=methodTitle)

        # -- inset axes to zoom
        axins = zoomed_inset_axes(ax, zoom=zoom, loc=4)
        axins.plot(self.data.raw['stress'], self.data.raw['work'], ls=':',
                   marker='o', lw=1.5, c='k', mfc='w')  # all data
        axins.plot(xRR, yRR, ls='-', c=colors[2], lw=1.125)
        axins.plot(sigmaRR, workRR, ls='', marker='+', c=colors[2])
        axins.plot(xCR, yCR, ls='-', c=colors[1], lw=1.125)
        if range2fitCR is not None or self.data.fitCc:
            axins.plot(sigmaCR, workCR, ls='', marker='x', c=colors[1])
        axins.plot(self.data.sigmaV, workSigmaV, ls='', marker='|', c='r',
                   ms=15, mfc='w', mew=1.5)
        axins.plot(self.sigmaP, self.wSigmaP, marker='o', c=colors[0], ms=7,
                   mfc='w', mew=1.5)
        # axins.spines['bottom'].set_visible(False)
        # axins.spines['right'].set_visible(False)
        axins.grid(False)
        axins.set(xlim=(-0.05*self.sigmaP, 1.25 * self.sigmaP),
                  ylim=(-0.05 * self.wSigmaP, 2.05 * self.wSigmaP),
                  xlabel=r'$\sigma^\prime_\mathrm{v}$ [kPa]',
                  ylabel=r'W [kJ m$^{-3}$]')
        axins.xaxis.tick_top()
        axins.xaxis.set_label_position('top')
        axins.yaxis.tick_right()
        axins.yaxis.set_label_position('right')
        axins.xaxis.set_tick_params(labelsize='small')
        axins.yaxis.set_tick_params(labelsize='small')
        axins.xaxis.set_minor_locator(mtick.AutoMinorLocator())
        axins.yaxis.set_minor_locator(mtick.AutoMinorLocator())
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
    >>> urlCSV = ''.join(['https://raw.githubusercontent.com/eamontoyaa/',
    >>>                   'data4testing/main/pysigmap/testData.csv'])
    >>> data = Data(pd.read_csv(urlCSV), sigmaV=75)
    >>> method = WangAndFrost(data)
    >>> method.getSigmaP(range2fitCR=None)
    >>> method.sigmaP, method.ocr
    (930.59421331855, 12.407922844247334)
    >>> method.getSigmaP(range2fitCR=[700, 10000])
    >>> method.sigmaP, method.ocr
    (718.9365936678746, 9.585821248904995)
    """

    def __init__(self, data):
        """Initialize the class."""
        # Heritage BeckerEtAl class
        BeckerEtAl.__init__(self, data)
        return

    def calculateDissipatedE(self, range2fitCR=None):
        """
        Calculate the accumulative dissipated strain energy corrected (ADSEC).

        Obtain the correction value based on the linear regression on the ADSE
        line or the tangent line to the steepest point of the cubic spline that
        passes through the points.

        Parameters
        ----------
        range2fitCR : list, tuple or array (length=2), optional
            Initial and final pressures between which the first order
            polynomial will be fit to the data on the compression range (CR)
            (range beyond the preconsolidation pressure). If None, the CR will
            be automatically fit to the same points used for calculating the
            compression index with the ``Data`` class only if it was calculated
            with a linear fit, otherwise, the steepest slope of the cubic
            spline that passes through the data will be used. The default is
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
        if range2fitCR is not None or self.data.fitCc:
            s2fitTE = self.data.cleaned['stress'][self.maskCR]
            w2fitTE = self.data.cleaned['ATSE'][self.maskCR]
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

    def getSigmaP(self, range2fitCR=None):
        """
        Return the value of the preconsolidation pressure or yield stress.

        Parameters
        ----------
        range2fitCR : list, tuple or array (length=2), optional
            Initial and final pressures between which the first order
            polynomial will be fit to the data on the compression range (CR)
            (range beyond the preconsolidation pressure). If None, the CR will
            be automatically fit to the same points used for calculating the
            compression index with the ``Data`` class only if it was calculated
            with a linear fit, otherwise, the steepest slope of the cubic
            spline that passes through the data will be used. The default is
            None.

        Returns
        -------
        fig : matplotlib figure
            Figure with the development of the method and the results.

        """
        # -- Post yield range or compression range (CR)
        self.maskCR = np.full(len(self.data.cleaned), False)
        if range2fitCR is not None or self.data.fitCc:  # Using a linear fit
            if range2fitCR is not None:
                idxInitCR = self.data.findStressIdx(
                    stress2find=range2fitCR[0], cleanedData=True)
                idxEndCR = self.data.findStressIdx(
                    stress2find=range2fitCR[1], cleanedData=True)
                self.maskCR[idxInitCR: idxEndCR] = True
            elif self.data.fitCc:
                self.maskCR = self.data.maskCc
            # -- Calculatting the dissipated strain energy corrected (ADSEC)
            self.calculateDissipatedE(range2fitCR)
            # -- Linear regresion of points on post yield line
            sigmaCR = self.data.cleaned['stress'][self.maskCR]
            disspE = self.data.cleaned['ADSEC'][self.maskCR]
            lcrInt, lcrSlope = polyfit(sigmaCR, disspE, deg=1)
            r2CR = r2_score(
                y_true=disspE, y_pred=polyval(sigmaCR, [lcrInt, lcrSlope]))

        else:  # Using the steepest point of a cubic spline
            # -- Calculatting the dissipated strain energy corrected (ADSEC)
            self.calculateDissipatedE(range2fitCR)
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

        xCR = np.linspace(self.sigmaP, self.data.cleaned['stress'].iloc[-1])
        yCR = polyval(xCR, [lcrInt, lcrSlope])

        # -- Plot compresibility curve
        fig = plt.figure(figsize=[9, 4.8])
        ax = fig.add_axes([0.08, 0.12, 0.55, 0.85])
        ax.plot(self.data.raw['stress'], self.data.raw['ATSEC'], ls=':',
                marker='v', lw=0.5, c='gray', mfc='w', label='Total energy')
        ax.plot(self.data.raw['stress'], self.data.raw['AESE'], ls=':',
                marker='s', lw=0.5, c='gray', mfc='w', label='Elastic energy')
        ax.plot(self.data.raw['stress'], self.data.raw['ADSEC'], ls=':',
                marker='o', lw=1.5, c='k', mfc='w', label='Dissipated energy')
        ax.hlines(y=abs(self.corrVal), xmin=0, xmax=self.sigmaP, lw=1.125,
                  color=colors[2])
        # Compression range
        ax.plot(xCR, yCR, ls='-', c=colors[1], lw=1.125,
                label='Compression range')
        if range2fitCR is not None or self.data.fitCc:
            ax.plot(sigmaCR, disspE, ls='', marker='x', c=colors[1],
                    label=f'Data for linear fit\n(R$^2={r2CR:.3f}$)')
        # Other plots
        ax.plot(self.data.sigmaV, -self.corrVal, ls='', marker='|', c='r',
                ms=15, mfc='w', mew=1.5, label=str().join([
                    r'$\sigma^\prime_\mathrm{v0}=$ ',
                    f'{self.data.sigmaV:.0f} kPa']))
        ax.plot(self.sigmaP, self.sseSigmaP, ls='', marker='o', c=colors[0],
                ms=7, mfc='w', mew=1.5,
                label=str().join([r'$\sigma^\prime_\mathrm{p}=$ ',
                                  f'{self.sigmaP:.0f} kPa\n',
                                  f'OCR= {self.ocr:.1f}']))
        # Other details
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set(ylabel='Specific strain energy [kPa]',
               xlabel=str().join(['Vertical effective stress, ',
                                  r'$\sigma^\prime_\mathrm{v}$ [kPa]']))
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.grid(False)
        ax.legend(bbox_to_anchor=(1.125, 0.5), loc=6,
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
