"""
data.py module.

Contains the class and its methods for determinig the preconsolidation
pressure from a consolidation test by the method proposed by
Casagrande (1963).
"""

# -- Required modules
import numpy as np
# from numpy.polynomial.polynomial import polyfit, polyval
# import matplotlib.pyplot as plt
import pandas as pd
# from sklearn.metrics import r2_score


# plt.rcParams['font.family'] = 'Serif'


class Data:
    """Data class."""

    def __init__(self, rawData, sigmaV, strainPercent=True, reloadStage=True,
                 doubleUnload=True):
        """
        Initialize the class.

        Parameters
        ----------
        rawData : pandas DataFrame
            Data from the consolidation test. It includes three series with
            the following and strict order: effective vertical stress,
            axial strain and void ratio.  The first row must include the
            initial void ratio of the sample.
        sigmaV : float
            Vertical effective stress of the tested sample.
        strainPercent : bool, optional
            Boolean to specify if the axial strain of the row data is in
            percent or not. The default is True.
        reloadStage : bool, optional
            Boolean to specify if the test included a reload stage. The default
            is True.
        doubleUnload : bool, optional
            Boolean to specify if the test included two reload stages. The
            default is True.

        Returns
        -------
        None.

        """
        self.raw = rawData
        self.sigmaV = sigmaV
        self.strainPercent = strainPercent
        self.reloadStage = reloadStage
        self.doubleUnload = doubleUnload
        self.preprocessing()
        self.getBreakIndices()
        self.clean()
        self.recompressionIdx()

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
            idx = -1
        else:
            data4finding = self.cleaned if cleanedData else self.raw
            idx = data4finding.query(f'stress >= {stress2find}').index[0]
        return idx

    def recompressionIdx(self):
        # -- Getting the recomression index
        deltaSigma = np.log10(self.raw['stress'][self.brkIdx2]) - \
            np.log10(self.raw['stress'][self.brkIdx1])
        delta_e = self.raw['e'][self.brkIdx2] - self.raw['e'][self.brkIdx1]
        self.idxCr = abs(delta_e / deltaSigma)
        return
