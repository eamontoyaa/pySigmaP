"""
``data.py`` module.

Contains the class and its methods for determinig the preconsolidation
pressure from a consolidation test by the method proposed by
Casagrande (1963).
"""

# -- Required modules
import numpy as np
from numpy.polynomial.polynomial import polyfit, polyval
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import r2_score
import matplotlib.ticker as mtick

plt.rcParams['font.family'] = 'Serif'
plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True


class Data:
    """
    ``Data`` class.

    When the object is instanced, the compression and recompression indices
    are automatically calculated with the default inputs of the methods.
    See the documentation of each method for more information regarding the
    default inputs.

    Attributes
    ----------
    rawData : pandas DataFrame
        Data from the consolidation test. It includes three series with the
        following and strict order: effective vertical stress, axial strain and
        void ratio.  The first row must include the initial void ratio of the
        sample.
    sigmaV : float
        Vertical effective stress of the tested sample.
    strainPercent : bool, optional
        Boolean to specify if the axial strain of the row data is in percent or
        not. The default is True.
    reloadStage : bool, optional
        Boolean to specify if the test included a reload stage. The default is
        True.
    doubleUnload : bool, optional
        Boolean to specify if the test included two reload stages. The default
        is True.

    Examples
    --------
    >>> data = Data(pd.read_csv('testData/testData.csv'), sigmaV=75)
    >>> # or
    >>> data = Data(pd.read_excel('testData/testData.xlsx'), sigmaV=75)
    >>> data.plot()
    >>> data.idxCc, data.idxCr
    (0.23829647651319996, 0.04873212611656529)
    >>> data.compressionIdx(range2fitCc=(1000, 8000))
    >>> data.plot()
    >>> data.idxCc
    0.22754962258271252
    >>> data.recompressionIdx(opt=2)
    >>> data.plot()
    >>> data.idxCr
    0.049481704199255676
    >>> data.recompressionIdx(opt=3)
    >>> data.plot()
    >>> data.idxCr
    0.05398842524225007
    """

    def __init__(self, rawData, sigmaV, strainPercent=True, reloadStage=True,
                 doubleUnload=True):
        """Initialize the class."""
        self.raw = rawData
        self.sigmaV = sigmaV
        self.strainPercent = strainPercent
        self.reloadStage = reloadStage
        self.doubleUnload = doubleUnload
        self.preprocessing()
        self.getBreakIndices()
        self.clean()
        self.compressionIdx()
        self.recompressionIdx()
        return

    def preprocessing(self):
        """
        Rename the series names and set the initial void ratio to an attribute.

        Returns
        -------
        None.
        """
        # Standardizing series names
        self.raw.columns = ['stress', 'strain', 'e']
        # Removing percentage format to strain values
        if self.strainPercent:
            self.raw['strain'] = self.raw['strain'].divide(100)
        # initial void ratio
        self.e_0 = self.raw['e'].iloc[0]
        return

    def getBreakIndices(self):
        """
        Find the break indices of the compressibility curve.

        The break indices are the following:

            - brkIdx1: The start of the first unload
            - brkIdx2: The end of the first unload (start of the first reload)
            - brkIdx3: Point on the NCL after first reload
            - brkIdx4: Index of the last point on the NCL

        Returns
        -------
        None.
        """
        for i in self.raw.index[:-1]:
            if self.raw['stress'][i+1] > self.raw['stress'][i] and \
                    self.raw['stress'][i+2] < self.raw['stress'][i+1]:
                brkIdx1 = i+1  # brkIdx1: start of the first unload
                break
        if self.reloadStage:
            for i in self.raw.index[brkIdx1+1:-1]:
                if self.raw['stress'][i+1] < self.raw['stress'][i] and \
                        self.raw['stress'][i+2] > self.raw['stress'][i+1]:
                    brkIdx2 = i+1  # brkIdx2: end of the first unload
                    break
            # brkIdx3: Point on the NCL after the first reload
            brkIdx3 = self.raw.query(f'stress == stress[{brkIdx1}]').index[1]
            # brkIdx4: index of the last point on the NCL
            brkIdx4 = self.raw.query('stress == stress.max()').index[0]
        else:
            brkIdx2 = self.raw.index[-1]
            brkIdx3 = None
            brkIdx4 = None

        self.brkIdx1 = brkIdx1
        self.brkIdx2 = brkIdx2
        self.brkIdx3 = brkIdx3
        self.brkIdx4 = brkIdx4
        return

    def clean(self):
        """
        Generate a cleaned copy of the raw data without unload-reload stages.

        Returns
        -------
        None.
        """
        if self.reloadStage:
            self.cleaned = pd.concat(
                [self.raw[0: self.brkIdx1+1],
                 self.raw[self.brkIdx3+1: self.brkIdx4+1]])
        else:
            self.cleaned = self.raw[0: self.brkIdx1+1]
        self.cleaned.reset_index(drop=True, inplace=True)  # update idx
        # -- Cubic spline that passes through the data
        sigmaLog = np.log10(self.cleaned['stress'][1:])
        cs = CubicSpline(x=sigmaLog, y=self.cleaned['e'][1:])
        self.eSigmaV = float(cs(np.log10(self.sigmaV)))  # void ratio at sigmaV
        return

    def findStressIdx(self, stress2find, cleanedData=True):
        """
        Return the ceiling index of a specified stress.

        Parameters
        ----------
        stress2find : float
            Stress whose ceiling index is wanted.
        cleanedData : bool, optional
            Boolean to indicate if the index must be find in the row or cleaned
            data. The default is True.

        Returns
        -------
        idx : int
            The ceiling index of the stress input.
        """
        if stress2find == 0:
            idx = 1
        elif stress2find > self.raw['stress'].max():
            idx = None
        else:
            data4finding = self.cleaned if cleanedData else self.raw
            idx = data4finding.query(f'stress >= {stress2find}').index[0]
        return idx

    def compressionIdx(self, range2fitCc=None):
        """
        Calculate the compression index (Cc).

        Parameters
        ----------
        range2fitCc : list, tuple or array (length=2), optional
            Initial and final pressures between which the first-order
            polynomial will be fit to the data on the normally consolidated
            line (NCL). If None, the Cc index will be calculated as the
            steepest slope of the cubic spline that passes through the data.
            The default is None.

        Returns
        -------
        None.
        """
        self.range2fitCc = range2fitCc
        if range2fitCc is None:  # Using the cubic spline for the NC line
            sigmaLog = np.log10(self.cleaned['stress'][1:])
            cs = CubicSpline(x=sigmaLog, y=self.cleaned['e'][1:])
            sigmaCS = np.linspace(sigmaLog.iloc[0], sigmaLog.iloc[-1], 500)
            steepestSlopeIdx = np.argmin(cs(sigmaCS, 1))
            idxCc = cs(sigmaCS, 1)[steepestSlopeIdx]
            idxCcInt = cs(sigmaCS[steepestSlopeIdx]) - \
                idxCc*sigmaCS[steepestSlopeIdx]
            self.fitCc = False
        else:  # Fitting a line to some points for the NC line
            idxInitCc = self.findStressIdx(
                stress2find=range2fitCc[0], cleanedData=True)
            idxEndCc = self.findStressIdx(
                stress2find=range2fitCc[1], cleanedData=True)
            maskCc = np.full(len(self.cleaned), False)
            maskCc[idxInitCc: idxEndCc] = True
            self.maskCc = maskCc
            # -- Linear regresion of points on normally consolidated line (Cc)
            sigmaCc = self.cleaned['stress'].iloc[maskCc]
            sigmaCclog = np.log10(sigmaCc)
            eCc = self.cleaned['e'].iloc[maskCc]
            idxCcInt, idxCc = polyfit(sigmaCclog, eCc, deg=1)
            r2Cc = r2_score(
                y_true=eCc, y_pred=polyval(sigmaCclog, [idxCcInt, idxCc]))
            self.r2Cc = r2Cc
            self.fitCc = True
        self.idxCc = abs(idxCc)
        self.idxCcInt = abs(idxCcInt)
        return

    def recompressionIdx(self, opt=1):
        """
        Calculate the recompression index (Cr).

        Parameters
        ----------
        opt : TYPE, optional
            Integer value to indicate which method will be used. Using only
            two points, the beginig and end of the unluad stage (opt=1); using
            all the points of the first unload stage (opt=2); using the points
            of the first unload and reload staged (opt=3). The default is 1.

        Returns
        -------
        None.
        """
        maskCr = np.full(len(self.raw), False)
        if opt == 1:
            maskCr[self.brkIdx1] = True
            maskCr[self.brkIdx2] = True
        elif opt == 2:
            maskCr[self.brkIdx1: self.brkIdx2+1] = True
        elif opt == 3:
            maskCr[self.brkIdx1: self.brkIdx3+1] = True
        # -- Linear regresion
        sigmaCr = self.raw['stress'].iloc[maskCr]
        sigmaCrlog = np.log10(sigmaCr)
        eCr = self.raw['e'].iloc[maskCr]
        idxCrInt, idxCr = polyfit(sigmaCrlog, eCr, deg=1)
        r2Cr = r2_score(
            y_true=eCr, y_pred=polyval(sigmaCrlog, [idxCrInt, idxCr]))
        self.maskCr = maskCr
        self.r2Cr = r2Cr
        self.idxCr = abs(idxCr)
        self.idxCrInt = idxCrInt
        return

    def plot(self):
        """
        Plot the compressibility curve, the Cc index and the Cr index.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Figure that includes the compressibility curve, the Cc index, the
            Cr index and the vertical effective stress of the sample tested.
        """
        # -- plotting
        fig = plt.figure(figsize=[9, 4.8])
        ax = fig.add_axes([0.08, 0.12, 0.55, 0.85])
        ax.plot(self.raw['stress'][1:], self.raw['e'][1:], ls='--', marker='o',
                lw=1, c='k', mfc='w', label='Experimental data')
        ax.plot(self.sigmaV, self.eSigmaV, ls='', marker='|', c='r', ms=15,
                mfc='w', label=str().join([r'$\sigma^\prime_\mathrm{v_0}=$ ',
                                           f'{self.sigmaV:.0f} kPa']))
        # Compression index
        x4Cc = np.linspace(
            self.cleaned['stress'].iloc[-4], self.cleaned['stress'].iloc[-1])
        y4Cc = -self.idxCc * np.log10(x4Cc) + self.idxCcInt
        ax.plot(x4Cc, y4Cc, ls='--', lw=0.8, color='darkgreen',
                label=str().join([r'$C_\mathrm{c}=$', f'{self.idxCc:.3f}']))
        if self.fitCc:
            ax.plot(self.cleaned['stress'].iloc[self.maskCc],
                    self.cleaned['e'].iloc[self.maskCc],
                    ls='', marker='x', lw=0.8, color='darkgreen',
                    label=f'Data for linear fit\n(R$^2={self.r2Cc:.3f}$)')
        # Recompression index
        x4Cr = np.linspace(self.raw['stress'].iloc[self.maskCr].min(),
                           self.raw['stress'].iloc[self.maskCr].max())
        y4Cr = -self.idxCr * np.log10(x4Cr) + self.idxCrInt
        ax.plot(x4Cr, y4Cr, ls='--', lw=0.8, color='darkred',
                label=str().join([r'$C_\mathrm{r}=$', f'{self.idxCr:.3f}']))
        ax.plot(self.raw['stress'].iloc[self.maskCr],
                self.raw['e'].iloc[self.maskCr], ls='', marker='+', lw=0.8,
                color='darkred',
                label=f'Data for linear fit\n(R$^2={self.r2Cr:.3f}$)')
        # other details
        ax.set(xscale='log', ylabel='Void ratio $(e)$',
               xlabel=str().join(['Vertical effective stress ',
                                  r'$(\sigma^\prime_\mathrm{v})$ [kPa]']))
        ax.xaxis.set_major_formatter(mtick.ScalarFormatter())
        ax.grid(True, which="both", ls='--', lw=0.5)
        ax.legend(bbox_to_anchor=(1.04, 0.5), loc=6,
                  title=r"\textbf{Compressibility curve}")
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
