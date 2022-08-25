"""
``data.py`` module.

Contains the class and its methods for loading and managing the data of the
compressibility curve.
"""

# -- Required modules
import numpy as np
from numpy.polynomial.polynomial import polyfit, polyval
from scipy.interpolate import CubicSpline
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas as pd
from mstools.mstools import r2_score

from pysigmap import figsize, colors


# plt.style.use("default")
mpl.rcParams.update(
    {
        "text.usetex": False,  # Use mathtext, not LaTeX
        "font.family": "serif",  # Use the Computer modern font
        "font.serif": "cmr10",
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "axes.unicode_minus": False,
    }
)


class Data:
    """``Data`` class.

    When the object is instanced, the compression and recompression indices
    are automatically calculated with the default inputs of the methods.
    See the documentation of each method for more information regarding the
    default inputs.

    Attributes
    ----------
    rawData : pandas DataFrame
        Data from the incremental loading oedometer test. It includes three
        series with the following and strict order: effective vertical stress,
        axial strain and void ratio. The first row must include the on-table
        void ratio of the specimen.
    sigmaV : float
        In situ effective vertical stress of the specimen.
    strainPercent : bool, optional
        Boolean to specify if the axial strain of the row data is in percent or
        not. The default is True.
    reloading : bool, optional
        Boolean to specify if the test included a reloading stage. The default
        is True.
    secondUnloading : bool, optional
        Boolean to specify if the test included two reloading stages. The
        default is True.

    Examples
    --------
    >>> urlCSV = ''.join(['https://raw.githubusercontent.com/eamontoyaa/',
    >>>                   'data4testing/main/pysigmap/testData.csv'])
    >>> urlXLSX = ''.join(['https://raw.githubusercontent.com/eamontoyaa/',
    >>>                    'data4testing/main/pysigmap/testData.xlsx'])
    >>> data = Data(pd.read_csv(urlCSV), sigmaV=75)
    >>> # or
    >>> data = Data(pd.read_excel(urlXLSX), sigmaV=75)
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

    def __init__(
        self,
        rawData,
        sigmaV,
        strainPercent=True,
        reloading=True,
        secondUnloading=True,
        use_this_UR_stage=0,
    ):
        """Initialize the class."""
        self.raw = rawData
        self.sigmaV = sigmaV
        self.strainPercent = strainPercent
        self.reloading = reloading
        self.secondUnloading = secondUnloading
        self.use_this_UR_stage = use_this_UR_stage
        self.preprocessing()
        self.getBreakIndices()
        self.clean()
        self.compressionIdx()
        self.recompressionIdx()
        return

    def preprocessing(self):
        """
        Rename the series names and create the on-table void ratio attribute.

        Returns
        -------
        None.
        """
        # Standardizing series names
        self.raw.columns = ["stress", "strain", "e"]
        # Removing percentage format to strain values
        if self.strainPercent:
            self.raw["strain"] = self.raw["strain"].divide(100)
        # On-table (initial) void ratio
        self.e_0 = self.raw["e"].iloc[0]
        return

    def getBreakIndices(self):
        """
        Find the break indices of the compressibility curve.

        The break indices are the following:

            - brkIdx1: Start of the first unloading stage
            - brkIdx2: End of the first unloading stage
            - brkIdx3: Point on the compression range after the first reloading
            - brkIdx4: Index of the last point on the compression range

        Returns
        -------
        None.
        """
        idx_unloading_init = [  # brkIdx1: Start of the first unloading
            i
            for i in self.raw.index[1:-1]
            if (
                self.raw["stress"][i] > self.raw["stress"][i - 1]
                and self.raw["stress"][i] > self.raw["stress"][i + 1]
            )
        ]
        if self.secondUnloading:
            idx_unloading_init.pop(-1)
        self.idx_unloading_init = idx_unloading_init
        brkIdx1 = idx_unloading_init[self.use_this_UR_stage]

        if self.reloading:  # Cases 2-3: Reloading and second unloading
            idx_reloading_init = [  # brkIdx2: End of the first unloading
                i
                for i in self.raw.index[1:-1]
                if (
                    self.raw["stress"][i] < self.raw["stress"][i - 1]
                    and self.raw["stress"][i] < self.raw["stress"][i + 1]
                )
            ]
            brkIdx2 = idx_reloading_init[self.use_this_UR_stage]
            self.idx_reloading_init = idx_reloading_init

            # brkIdx3: Points on the NCL after the each reloading
            idx_reloading_ncl = [
                self.raw.query(f"stress == stress[{i}]").index[1]
                for i in idx_unloading_init
            ]
            brkIdx3 = idx_reloading_ncl[self.use_this_UR_stage]
            self.idx_reloading_ncl = idx_reloading_ncl
        else:  # Case 1: No reloading - No second unloading
            brkIdx2 = self.raw.index[-1]  # last point of the unluading
            brkIdx3 = None
        # brkIdx4: index of the last point on the NCL
        brkIdx4 = self.raw.query("stress == stress.max()").index[0]

        self.brkIdx1 = brkIdx1
        self.brkIdx2 = brkIdx2
        self.brkIdx3 = brkIdx3
        self.brkIdx4 = brkIdx4
        return

    def clean(self):
        """
        Clean the raw data removing the unloading and reloading stages.

        Returns
        -------
        None.
        """
        if self.reloading:
            idx_init = [-1] + self.idx_reloading_ncl
            idx_end = self.idx_unloading_init + [self.brkIdx4]
            self.cleaned = pd.concat(
                [
                    self.raw[idx_init[i] + 1 : idx_end[i] + 1]
                    for i in range(len(idx_init))
                ]
            )
        else:
            self.cleaned = self.raw[: self.brkIdx1 + 1]
        self.cleaned.reset_index(drop=True, inplace=True)
        sigmaLog = np.log10(self.cleaned["stress"][1:])
        cs = CubicSpline(x=sigmaLog, y=self.cleaned["e"][1:])
        self.eSigmaV = float(cs(np.log10(self.sigmaV)))
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
            return 1
        elif stress2find > self.raw["stress"].max():
            return None
        else:
            data4finding = self.cleaned if cleanedData else self.raw
            return data4finding.query(f"stress >= {stress2find}").index[0]

    def compressionIdx(self, range2fitCc=None):
        """
        Calculate the compression index (Cc).

        Parameters
        ----------
        range2fitCc : list, tuple or array (length=2), optional
            Initial and final pressures between which the first-order
            polynomial will be fit to the data on the compression range. If
            None, the Cc index will be calculated as the
            steepest slope of the cubic spline that passes through the data.
            The default is None.

        Returns
        -------
        None.
        """
        self.range2fitCc = range2fitCc
        if range2fitCc is None:  # Using a cubic spline
            sigmaLog = np.log10(self.cleaned["stress"][1:])
            cs = CubicSpline(x=sigmaLog, y=self.cleaned["e"][1:])
            sigmaCS = np.linspace(sigmaLog.iloc[0], sigmaLog.iloc[-1], 500)
            steepestSlopeIdx = np.argmin(cs(sigmaCS, 1))
            idxCc = cs(sigmaCS, 1)[steepestSlopeIdx]
            idxCcInt = (
                cs(sigmaCS[steepestSlopeIdx])
                - idxCc * sigmaCS[steepestSlopeIdx]
            )
            self.fitCc = False
        else:  # Using a linear fit of a first-order polynomial
            idxInitCc = self.findStressIdx(
                stress2find=range2fitCc[0], cleanedData=True
            )
            idxEndCc = self.findStressIdx(
                stress2find=range2fitCc[1], cleanedData=True
            )
            maskCc = np.full(len(self.cleaned), False)
            maskCc[idxInitCc:idxEndCc] = True
            self.maskCc = maskCc
            # -- Linear regresion
            sigmaCc = self.cleaned["stress"].iloc[maskCc]
            sigmaCclog = np.log10(sigmaCc)
            eCc = self.cleaned["e"].iloc[maskCc]
            idxCcInt, idxCc = polyfit(sigmaCclog, eCc, deg=1)
            r2Cc = r2_score(
                y_true=eCc, y_pred=polyval(sigmaCclog, [idxCcInt, idxCc])
            )
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
            two points, the start and end of the unluading stage (opt=1); using
            all the points of the first unloading stage (opt=2); using the
            points of the first unloading and reloading stages (opt=3). The
            default is 1.

        Returns
        -------
        None.
        """
        maskCr = np.full(len(self.raw), False)
        if opt == 1:
            maskCr[self.brkIdx1] = True
            maskCr[self.brkIdx2] = True
        elif opt == 2:
            maskCr[self.brkIdx1 : self.brkIdx2 + 1] = True
        elif opt == 3:
            maskCr[self.brkIdx1 : self.brkIdx3 + 1] = True
        # -- Linear regresion
        sigmaCr = self.raw["stress"].iloc[maskCr]
        sigmaCrlog = np.log10(sigmaCr)
        eCr = self.raw["e"].iloc[maskCr]
        idxCrInt, idxCr = polyfit(sigmaCrlog, eCr, deg=1)
        r2Cr = r2_score(
            y_true=eCr, y_pred=polyval(sigmaCrlog, [idxCrInt, idxCr])
        )
        self.maskCr = maskCr
        self.r2Cr = r2Cr
        self.idxCr = abs(idxCr)
        self.idxCrInt = idxCrInt
        return

    def plot(self):
        """
        Plot the compressibility curve, Cc and Cr indices.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Figure that includes the compressibility curve, Cc, Cr and the in
            situ effective vertical stress of the specimen.
        """
        # -- plotting
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.08, 0.12, 0.55, 0.85])
        ax.plot(
            self.raw["stress"][1:],
            self.raw["e"][1:],
            ls=(0, (1, 1)),
            marker="o",
            lw=1.5,
            c="k",
            mfc="w",
            label="Experimental data",
        )
        ax.plot(
            self.sigmaV,
            self.eSigmaV,
            ls="",
            marker="|",
            c="r",
            ms=15,
            mfc="w",
            mew=1.5,
            label=str().join(
                ["$\sigma^\prime_\mathrm{v_0}=$ ", f"{self.sigmaV:.0f} kPa"]
            ),
        )
        # Compression index
        x4Cc = np.linspace(
            self.cleaned["stress"].iloc[-4], self.cleaned["stress"].iloc[-1]
        )
        y4Cc = -self.idxCc * np.log10(x4Cc) + self.idxCcInt
        ax.plot(
            x4Cc,
            y4Cc,
            ls="-",
            lw=1.125,
            color=colors[1],
            label=str().join(["$C_\mathrm{c}=$", f"{self.idxCc:.3f}"]),
        )
        if self.fitCc:
            ax.plot(
                self.cleaned["stress"].iloc[self.maskCc],
                self.cleaned["e"].iloc[self.maskCc],
                ls="",
                marker="x",
                color=colors[1],
                label=f"Data for linear fit\n(R$^2={self.r2Cc:.3f}$)",
            )
        # Recompression index
        x4Cr = np.linspace(
            self.raw["stress"].iloc[self.maskCr].min(),
            self.raw["stress"].iloc[self.maskCr].max(),
        )
        y4Cr = -self.idxCr * np.log10(x4Cr) + self.idxCrInt
        ax.plot(
            x4Cr,
            y4Cr,
            ls="-",
            lw=1.125,
            color=colors[2],
            label=str().join(["$C_\mathrm{r}=$", f"{self.idxCr:.3f}"]),
        )
        ax.plot(
            self.raw["stress"].iloc[self.maskCr],
            self.raw["e"].iloc[self.maskCr],
            ls="",
            marker="+",
            color=colors[2],
            label=f"Data for linear fit\n(R$^2={self.r2Cr:.3f}$)",
        )
        # other details
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set(
            xscale="log",
            ylabel="Void ratio, $e$",
            xlabel=str().join(
                [
                    "Effective vertical stress, ",
                    "$\sigma^\prime_\mathrm{v}$ [kPa]",
                ]
            ),
        )
        ax.xaxis.set_major_formatter(mtick.ScalarFormatter())
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.grid(False)
        ax.legend(
            bbox_to_anchor=(1.125, 0.5),
            loc=6,
            title="$\\bf{Compressibility\ curve}$",
        )
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
