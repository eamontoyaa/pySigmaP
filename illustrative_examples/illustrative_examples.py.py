"""
appendix.py module.

This script contains the code used to generate the figures related to the
section 'Software Functionalities and Illustrative Examples' of the paper
'An Open Source Application Software to Determine the
Preconsolidation Pressure of Soils in Incremental Loading Oedometer Testing:
pySigmaP' by Exneyder A. Montoya-Araque, A. J. Aparicio-Ortube, David G.
Zapata-Medina and  Luis G. Arboleda-Monsalve.
"""

# Block 1: Input loading data from an external file
from pysigmap.energy import WangAndFrost
from pysigmap.energy import BeckerEtAl
from pysigmap.bilog import Bilog
from pysigmap.boone import Boone
from pysigmap.pachecosilva import PachecoSilva
from pysigmap.casagrande import Casagrande
import pandas as pd
from pysigmap.data import Data
url = ''.join([
    'https://raw.githubusercontent.com/',
    'eamontoyaa/data4testing/',
    'main/pysigmap/testData.csv'])
data = Data(pd.read_csv(url), sigmaV=75)
fig = data.plot()  # Figure 2a

# Block 2: Cc and Cr calculated following published criteria
# 2.1 - Default parameters: Cc (maximum slope) – Cr (opt=1)
data.compressionIdx(range2fitCc=None)
data.recompressionIdx(opt=1)
fig = data.plot()  # Figure 2a
# 2.2: Cc (two last points) – Cr (opt=2)
data.compressionIdx(range2fitCc=(3000, 8000))
data.recompressionIdx(opt=2)
fig = data.plot()  # Figure 2b
# 2.3: Cc (four last points) – Cr (opt=3)
data.compressionIdx(range2fitCc=(700, 8000))
data.recompressionIdx(opt=3)
fig = data.plot()  # Figure 2c

# Block 3: Computation of 〖σ'〗_"p"  via the Casagrande method
method = Casagrande(data)
# 3.1 - default parameters: cubic spline function
fig = method.getSigmaP(mcp=None, range2fitFOP=None, loglog=True)  # Figure 3a
# 3.2: fourth order polynomial (FOP)
fig = method.getSigmaP(range2fitFOP=[20, 5000], loglog=True)  # Figure 3b
# 3.3: MCP manually introduced
fig = method.getSigmaP(mcp=200)  # Not shown

# Block 4: Computation of 〖σ'〗_"p"  via the Pacheco Silva and Boone methods
# 4.1: Pacheco Silva method
method = PachecoSilva(data)
fig = method.getSigmaP()  # Figure 3c
# 4.2: Boone method
method = Boone(data)
fig = method.getSigmaP()  # Figure 3d

# Block 5: Computation of 〖σ'〗_"p"  via the bilogarithmic methods
method = Bilog(data)
# 5.1: Butterfield method
fig = method.getSigmaP(range2fitRR=None, range2fitCR=None, opt=1)  # Figure 4a

# 5.2: Oikawa method
fig = method.getSigmaP(
    range2fitRR=None, range2fitCR=[1000, 5000], opt=2)  # Figure 4b
# 5.3: Onitsuka et al. method
fig = method.getSigmaP(
    range2fitRR=[0, 30], range2fitCR=[1000, 9000], opt=3)  # Figure 4c

# Block 6: Computation of 〖σ'〗_"p"  via the strain energy methods
method = BeckerEtAl(data)
# 6.1: Becker et al. method
fig = method.getSigmaP(range2fitRR=None, range2fitCR=None,
                       morinFormulation=False, zoom=5.5)  # Figure 5a
# 6.2: Morin method
fig = method.getSigmaP(range2fitRR=[0, 100], range2fitCR=[700, 9000],
                       morinFormulation=True, zoom=5.5)  # Figure 5b
# 6.3: Wang and Frost method
method = WangAndFrost(data)
fig = method.getSigmaP(range2fitCR=None)  # Figure 5c
