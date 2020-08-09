"""
beckeretal.py module.

Contains the class and its methods for determinig the preconsolidation
pressure from a consolidation test by the energy methods proposed by
Becker et al. (1987) and Wang & Frost (2004).

References
----------
Becker, D. E., Crooks, J. H. A., Been, K., & Jefferies, M. G. (1987). Work as a
criterion for determining in situ and yield stresses in clays. Canadian
Geotechnical Journal, 24, 4, 549-564, https://doi.org/10.1139/t87-070

Wang, L. B., & Frost, J. D. (2004). Dissipated strain energy method for
determining preconsolidation pressure. Canadian Geotechnical Journal, 41, 4,
760-768, https://doi.org/10.1139/t04-013

"""

# -- Required modules
import numpy as np
from numpy.polynomial.polynomial import polyfit, polyval
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from sklearn.metrics import r2_score

plt.rcParams['font.family'] = 'Serif'
plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True


class BeckerEtAl():
    """BeckerEtAl class."""

    def __init__(self, data, morinApproach=False):
        """
        Initialize the BeckerEtAl class.

        Instance an object to perform the energy method by Becker et al. (1987)
        for determining the preconsolidation pressure from a unidimensional
        consolidation test.

        Parameters
        ----------
        data : Object instanced from the Data class.
            Contains the data structure from the consolidation test.

        Returns
        -------
        None.

        """
        self.data = data
        self.morinApproach = morinApproach
        # Applying intinsec method
        self.calculateWork()
        return

    def calculateWork(self):
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
        if self.morinApproach:  # Work per unit volume of solids
            deltaWork *= (1 + self.data.e_0)
        self.data.raw['deltaWork'] = deltaWork
        self.data.raw['work'] = np.cumsum(deltaWork)
        self.data.clean()  # Data without unloads
        return

    def getSigmaP(self, range2fitNCL=None, zoom=2):
        """
        Return the value of the preconsolidation pressure or yield stress.

        Parameters
        ----------
        range2fitNCL : list, tuple or array (length=2), optional
            Initial and final pressures between which the first-order
            polynomial is fit to the compressibility curve on the normally
            consolidated line (NCL). If None, the NCL is fit to the last three
            points of the curve. The default is None.
        zoom : int, optional
            Zoom value to get close to the firsts points of the test and the
            preconsolidation pressure. The default is -3.

        Returns
        -------
        fig : matplotlib figure
            Figure with the development of the method and the results.

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

        # -- Linear regresion of points on the preyield line (PYL)
        sigmaPYL = self.data.cleaned['stress'][0: sigmaVidx+1]
        workPYL = self.data.cleaned['work'][0: sigmaVidx+1]
        p1_0, p1_1 = polyfit(sigmaPYL, workPYL, deg=1)
        r2PYL = r2_score(
            y_true=workPYL, y_pred=polyval(sigmaPYL, [p1_0, p1_1]))
        xPYL = np.linspace(0, self.data.cleaned['stress'].iloc[idxInitNCL])
        yPYL = polyval(xPYL, [p1_0, p1_1])

        # -- Linear regresion of points on the normally consolidated line (NCL)
        sigmaNCL = self.data.cleaned['stress'][idxInitNCL: idxEndNCL]
        workNCL = self.data.cleaned['work'][idxInitNCL: idxEndNCL]
        p2_0, p2_1 = polyfit(sigmaNCL, workNCL, deg=1)
        r2NCL = r2_score(
            y_true=workNCL, y_pred=polyval(sigmaNCL, [p2_0, p2_1]))
        xNCL = np.linspace(-p2_0/p2_1, self.data.cleaned['stress'].iloc[-1])
        yNCL = polyval(xNCL, [p2_0, p2_1])

        # -- Preconsolitadion pressure
        workSigmaV = polyval(self.data.sigmaV, [p1_0, p1_1])
        self.sigmaP = (p2_0 - p1_0) / (p1_1 - p2_1)
        w_sigmaP = polyval(self.sigmaP, [p1_0, p1_1])

        # -- plot compresibility curve
        fig = plt.figure(figsize=[9, 4.8])
        # ax = fig.add_subplot(111)
        ax = fig.add_axes([0.08, 0.12, 0.65, 0.85])
        ax.plot(self.data.raw['stress'], self.data.raw['work'], ls='--',
                marker='o', lw=1, c='k', mfc='w', label='Data')  # all data
        ax.plot(xPYL, yPYL, ls='--', c='darkcyan', lw=0.8,
                label=f'Preyield line \n(R$^2={r2PYL:.3f}$)')
        ax.plot(sigmaPYL, workPYL, ls='', marker='.', c='darkcyan')
        ax.plot(xNCL, yNCL, ls='--', c='crimson', lw=0.8,
                label=f'Postyield line \n(R$^2={r2NCL:.3f}$)')
        ax.plot(sigmaNCL, workNCL, ls='', marker='.', c='crimson')
        ax.plot(self.data.sigmaV, workSigmaV, ls='', marker='|', c='r', ms=15,
                mfc='w', label=str().join([r'$\sigma^\prime_\mathrm{v0}=$ ',
                                           f'{self.data.sigmaV:.0f} kPa']))
        ax.plot(self.sigmaP, w_sigmaP, ls='', marker='D', c='r', ms=5, mfc='w',
                label=str().join([r'$\sigma^\prime_\mathrm{p}=$ ',
                                  f'{self.sigmaP:.0f} kPa']))
        # other details
        ax.set(ylabel='Total work per unit vol. (W) [kJ m$^{-3}$]',
               xlabel=str().join(['Vertical effective stress ',
                                  r'$(\sigma^\prime_\mathrm{v})$ [kPa]']))
        ax.grid(True, ls='--', lw=0.5)
        methodTitle = r"\textbf{Becker \textit{et al.}'s method}"
        if self.morinApproach:
            methodTitle = r"\textbf{Morin's method}"
        ax.legend(bbox_to_anchor=(1.02, 0.5), loc=6, title=methodTitle)

        # -- inset axes to zoom
        axins = zoomed_inset_axes(ax, zoom=zoom, loc=4)
        axins.plot(self.data.raw['stress'], self.data.raw['work'], ls='--',
                   marker='o', lw=1, c='k', mfc='w')  # all data
        axins.plot(xPYL, yPYL, ls='--', c='darkcyan', lw=0.8)
        axins.plot(sigmaPYL, workPYL, ls='', marker='.', c='darkcyan')
        axins.plot(xNCL, yNCL, ls='--', c='crimson', lw=0.8)
        axins.plot(sigmaNCL, workNCL, ls='', marker='.', c='crimson')
        axins.plot(self.data.sigmaV, workSigmaV, ls='', marker='|', c='r',
                   ms=15, mfc='w')
        axins.plot(self.sigmaP, w_sigmaP, marker='D', c='r', ms=5, mfc='w')
        axins.grid(True, ls='--', lw=0.5)
        axins.set(xlim=(-0.05*self.sigmaP, 1.25 * self.sigmaP),
                  ylim=(-0.05 * w_sigmaP, 2.05 * w_sigmaP),
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
    """WangAndFrost class."""

    def __init__(self, data):
        """
        Initialize the WangAndFrost class.

        Instance an object to perform the energy method by Becker et al. (1987)
        for determining the preconsolidation pressure from a unidimensional
        consolidation test.

        Parameters
        ----------
        data : Object instanced from the Data class.
            Contains the data structure from the consolidation test.

        Returns
        -------
        None.

        """
        # initializing BeckerEtAl class
        BeckerEtAl.__init__(self, data)
        return

    def calculateDissipatedE(self, idxInitNCL=-3, idxEndNCL=None):
        """
        Calculate the accumulative dissipated strain energy corrected.

        Obtain the correction value based on the  linear regression on the ADSE
            line.

        Parameters
        ----------
        idxInitNCL : int, optional
            Initial index of the cleaned data for the linear regression on the
            ADSE line. The default is -3.
        idxEndNCL : int, optional
            final index of the cleaned data for the linear regression on the
            ADSE line. The default is None.

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
        s2fitTE = self.data.cleaned['stress'][idxInitNCL:idxEndNCL]
        w2fitTE = self.data.cleaned['ATSE'][idxInitNCL:idxEndNCL]
        self.corrVal, _ = polyfit(s2fitTE, w2fitTE, deg=1)
        # -- Calculating the accumulative total strain energy corrected (ATSEC)
        self.data.raw['ATSEC'] = self.data.raw['ATSE'] + abs(self.corrVal)
        # -- Accumulative dissipated strain energy corrected (ADSEC)
        self.data.raw['ADSEC'] = self.data.raw['ADSE'] + abs(self.corrVal)
        # -- Updating data without unloads
        self.data.clean()
        return

    def getSigmaP(self, range2fitNCL=None):
        """
        Return the value of the preconsolidation pressure or yield stress.

        Parameters
        ----------
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
        if range2fitNCL is None:  # Indices for fitting the NCL line
            idxInitNCL = -3
            idxEndNCL = None
        else:
            idxInitNCL = self.data.findStressIdx(
                stress2find=range2fitNCL[0], cleanedData=True)
            idxEndNCL = self.data.findStressIdx(
                stress2find=range2fitNCL[1], cleanedData=True) if \
                range2fitNCL[1] < self.data.cleaned['stress'].max() else None
        # -- Calculatting the dissipated strain energy corrected (ADSEC)
        self.calculateDissipatedE(idxInitNCL, idxEndNCL)

        # -- Linear regresion of points on the ADSEC line
        sigmaADSEC = self.data.cleaned['stress'][idxInitNCL:idxEndNCL]
        disspE = self.data.cleaned['ADSEC'][idxInitNCL:idxEndNCL]
        p1_0, p1_1 = polyfit(sigmaADSEC, disspE, deg=1)
        r2ADSEC = r2_score(
            y_true=disspE, y_pred=polyval(sigmaADSEC, [p1_0, p1_1]))
        xADSEC = np.linspace(0, self.data.cleaned['stress'].iloc[-1])
        yADSEC = polyval(xADSEC, [p1_0, p1_1])

        # Preconsolitadion pressure
        self.sigmaP = (abs(self.corrVal) - p1_0) / p1_1
        w_sigmaP = polyval(self.sigmaP, [p1_0, p1_1])

        # -- plot compresibility curve
        fig = plt.figure(figsize=[9, 4.8])
        # ax = fig.add_subplot(111)
        ax = fig.add_axes([0.08, 0.12, 0.65, 0.85])
        ax.plot(self.data.raw['stress'], self.data.raw['ATSEC'], ls=':',
                marker='v', lw=0.3, c='gray', mfc='w', label='Total energy')
        ax.plot(self.data.raw['stress'], self.data.raw['AESE'], ls=':',
                marker='s', lw=0.3, c='gray', mfc='w', label='Elastic energy')
        ax.plot(self.data.raw['stress'], self.data.raw['ADSEC'], ls='--',
                marker='o', lw=1, c='k', mfc='w', label='Dissipated energy')
        ax.hlines(y=abs(self.corrVal), xmin=0, xmax=self.sigmaP, lw=0.8,
                  color='darkcyan')
        ax.plot(xADSEC, yADSEC, ls='--', c='crimson', lw=0.8,
                label=f'Linear regression \n(R$^2={r2ADSEC:.3f}$)')
        ax.plot(sigmaADSEC, disspE, ls='', marker='.', c='crimson')
        ax.plot(self.sigmaP, w_sigmaP, ls='', marker='D', c='r', ms=5,
                mfc='w', label=str().join([r'$\sigma^\prime_\mathrm{p}=$ ',
                                           f'{self.sigmaP:.0f} kPa']))
        # other details
        ax.set(ylabel='Specific strain energy [kJ m$^{-3}$]',
               xlabel=str().join(['Vertical effective stress ',
                                  r'$(\sigma^\prime_\mathrm{v})$ [kPa]']))
        ax.grid(True, ls='--', lw=0.5)
        ax.legend(bbox_to_anchor=(1.02, 0.5), loc=6,
                  title=r"\textbf{Wang \& Frost's method}")
        return fig


# %%
"""
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
"""
