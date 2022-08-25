"""
``bilog.py`` module
-------------------

Contains the class and its methods to determine the preconsolidation
pressure from a compressibility curve via the Bilogarithmic methods proposed
by Butterfield (1979), Oikawa (1987) and Onitsuka et al. (1995).

References
~~~~~~~~~~
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
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from mstools.mstools import r2_score

from pysigmap import figsize, colors

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


class Bilog:
    """
    ``Bilog`` class.

    When the object is instanced, the method ``getSigmaP()`` calculates the
    preconsolidation pressure by the methods proposed by Butterfield (1979),
    Oikawa (1987) and Onitsuka et al. (1995) based on the parameters of the
    method. See the method documentation for more information.

    Attributes
    ----------
    data : Object instanced from the ``Data`` class
        Contains the data structure from the oedometer test. See the class
        documentation for more information.

    Examples
    --------
    >>> urlCSV = ''.join(['https://raw.githubusercontent.com/eamontoyaa/',
    >>>                   'data4testing/main/pysigmap/testData.csv'])
    >>> data = Data(pd.read_csv(urlCSV), sigmaV=75)
    >>> method = Bilog(data)
    >>> method.getSigmaP(opt=1)  # Butterfield (1979)
    >>> method.sigmaP, method.ocr
    (433.29446692547555, 5.777259559006341)
    >>> method.getSigmaP(range2fitCR=[1000, 9000], opt=1) # Butterfield
    >>> method.sigmaP, method.ocr
    (399.3361458736937, 5.324481944982582)
    >>> method.getSigmaP(opt=2)  # Oikawa (1987)
    >>> method.sigmaP, method.ocr
    (433.2944669254731, 5.777259559006308)
    >>> method.getSigmaP(range2fitCR=[1000, 9000], opt=2)  # Oikawa
    >>> method.sigmaP, method.ocr
    (399.3361458736935, 5.32448194498258)
    >>> method.getSigmaP(opt=3)  # Onitsuka et al. (1995)
    >>> method.sigmaP, method.ocr
    (433.2944669254753, 5.777259559006338)
    >>> method.getSigmaP(range2fitCR=[1000, 9000], opt=3)  # Onitsuka et al.
    >>> method.sigmaP, method.ocr
    (399.3361458736939, 5.324481944982585)
    """

    def __init__(self, data):
        """Initialize the class."""
        self.data = data
        return

    def getSigmaP(self, range2fitRR=None, range2fitCR=None, opt=1):
        """
        Return the value of the preconsolidation pressure.

        Parameters
        ----------
        range2fitRR : list, tuple or array (length=2), optional
            Initial and final pressures between which the first order
            polynomial will be fitted to the data on the recompression range
            (RR) (range before the preconsolidation pressure). If None, the
            first order polynomial will be fitted from the first point of the
            curve to the point before the in situ effective vertical stress.
            The default is None.
        range2fitCR : list, tuple or array (length=2), optional
            Initial and final pressures between which the first order
            polynomial will be fitted to the data on the compression range (CR)
            (range above the preconsolidation pressure). If None, the CR will
            be automatically fitted to the same points used for calculating the
            compression index with the ``Data`` class only if it was calculated
            with a linear fit, otherwise, the steepest slope of the cubic
            spline that passes through the data will be used. The default is
            None.
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

        self.data.raw["vol"] = self.data.raw["e"] + 1
        self.data.clean()  # -- Updating data without unloads

        # def ticks(x, pos): return f'$e^{np.log(x):.0f}$'

        if range2fitRR is None:  # Indices for fitting the RR line
            idxInitRR = 1
            idxEndRR = (
                self.data.findStressIdx(
                    stress2find=self.data.sigmaV, cleanedData=True
                )
                - 1
            )
        else:
            idxInitRR = self.data.findStressIdx(
                stress2find=range2fitRR[0], cleanedData=True
            )
            idxEndRR = (
                self.data.findStressIdx(
                    stress2find=range2fitRR[1], cleanedData=True
                )
                - 1
            )

        # -- Linear regresion of points on the recompression range(RR)
        sigmaRR = self.data.cleaned["stress"][idxInitRR : idxEndRR + 1]
        volRR = self.data.cleaned["vol"][idxInitRR : idxEndRR + 1]
        sigmaRRlog = transformX(sigmaRR, opt)
        volRRlog = transformY(volRR, opt)
        p1_0, p1_1 = polyfit(sigmaRRlog, volRRlog, deg=1)
        r2RR = r2_score(
            y_true=volRRlog, y_pred=polyval(sigmaRRlog, [p1_0, p1_1])
        )

        # -- Cubic spline that passes through the data
        sigmaLog = transformX(self.data.cleaned["stress"][1:], opt)
        volLog = transformY(self.data.cleaned["vol"][1:], opt)
        cs = CubicSpline(x=sigmaLog, y=volLog)
        # Specific volume at sigma V
        volSigmaV = cs(transformX(self.data.sigmaV, opt))

        # -- Compression range (CR)
        self.maskCR = np.full(len(self.data.cleaned), False)
        if range2fitCR is not None or self.data.fitCc:  # Using a linear fit
            if range2fitCR is not None:
                idxInitCR = self.data.findStressIdx(
                    stress2find=range2fitCR[0], cleanedData=True
                )
                idxEndCR = self.data.findStressIdx(
                    stress2find=range2fitCR[1], cleanedData=True
                )
                self.maskCR[idxInitCR:idxEndCR] = True
            elif self.data.fitCc:
                self.maskCR = self.data.maskCc
            # -- Linear regresion of points on post yield line
            sigmaCR = self.data.cleaned["stress"][self.maskCR]
            sigmaCRlog = transformX(sigmaCR, opt)
            volCR = self.data.cleaned["vol"][self.maskCR]
            volCRlog = transformY(volCR, opt)
            lcrInt, lcrSlope = polyfit(sigmaCRlog, volCRlog, deg=1)
            r2CR = r2_score(
                y_true=volCRlog, y_pred=polyval(sigmaCRlog, [lcrInt, lcrSlope])
            )

        else:  # Using the steepest point of a cubic spline
            sigmaCS = np.linspace(sigmaLog.iloc[0], sigmaLog.iloc[-1], 500)
            steepestSlopeIdx = np.argmin(cs(sigmaCS, 1))
            lcrSlope = cs(sigmaCS, 1)[steepestSlopeIdx]
            lcrInt = (
                cs(sigmaCS[steepestSlopeIdx])
                - lcrSlope * sigmaCS[steepestSlopeIdx]
            )

        xScl = np.e if opt == 1 else 10  # log scale for plotting
        yLabel = "$\log_{10} (1+e)$" if opt == 2 else "$\ln (1+e)$"

        # -- Preconsolitadion pressure
        self.sigmaP = xScl ** ((lcrInt - p1_0) / (p1_1 - lcrSlope))
        self.vSigmaP = polyval(transformX(self.sigmaP, opt), [p1_0, p1_1])
        self.ocr = self.sigmaP / self.data.sigmaV

        # -- Lines of the bilogarithmic methods
        xCR = np.linspace(self.sigmaP, self.data.cleaned["stress"].iloc[-1])
        yCR = polyval(transformX(xCR, opt), [lcrInt, lcrSlope])

        # Recompression range
        xRRFit = np.linspace(sigmaRR.iloc[0], self.sigmaP)
        yRRFit = polyval(transformX(xRRFit, opt), [p1_0, p1_1])

        # -- plot compresibility curve
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.08, 0.12, 0.55, 0.85])
        ax.semilogx(
            self.data.raw["stress"][1:],
            transformY(self.data.raw["vol"][1:], opt),
            base=xScl,
            ls=(0, (1, 1)),
            marker="o",
            lw=1.5,
            c="k",
            mfc="w",
            label="Compressibility curve",
        )
        methods = ["Butterfield", "Oikawa", "Onitsuka\ et\ al."]
        # Recompression range
        ax.plot(
            xRRFit,
            yRRFit,
            ls="-",
            c=colors[2],
            lw=1.125,
            label="Recompression range",
        )
        ax.plot(
            sigmaRR,
            volRRlog,
            ls="",
            marker="+",
            c=colors[2],
            label=f"Data for linear fit\n(R$^2={r2RR:.3f}$)",
        )
        # Compression range
        ax.plot(
            xCR, yCR, ls="-", c=colors[1], lw=1.125, label="Compression range"
        )
        if range2fitCR is not None or self.data.fitCc:
            ax.plot(
                sigmaCR,
                volCRlog,
                ls="",
                marker="x",
                c=colors[1],
                label=f"Data for linear fit\n(R$^2={r2CR:.3f}$)",
            )
        # Other plots
        ax.plot(
            self.data.sigmaV,
            volSigmaV,
            ls="",
            marker="|",
            c="r",
            ms=15,
            mfc="w",
            mew=1.5,
            label=str().join(
                [
                    "$\sigma^\prime_\mathrm{v0}=$ ",
                    f"{self.data.sigmaV:.0f} kPa",
                ]
            ),
        )
        ax.plot(
            self.sigmaP,
            self.vSigmaP,
            ls="",
            marker="o",
            c=colors[0],
            ms=7,
            mfc="w",
            mew=1.5,
            label=str().join(
                [
                    "$\sigma^\prime_\mathrm{p}=$ ",
                    f"{self.sigmaP:.0f} kPa\n",
                    f"OCR= {self.ocr:.1f}",
                ]
            ),
        )
        # Other details
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set(
            ylabel=f"Specific volume, {yLabel}",
            xlabel=str().join(
                [
                    "Effective vertical stress, ",
                    "$\sigma^\prime_\mathrm{v}$ [kPa]",
                ]
            ),
        )
        ax.xaxis.set_major_formatter(mtick.ScalarFormatter())
        # ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.grid(False)
        ax.legend(
            bbox_to_anchor=(1.125, 0.5),
            loc=6,
            title=str().join(["$\\bf{", f"{methods[opt-1]}", "\ method}$"]),
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
