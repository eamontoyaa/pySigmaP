"""test.py module."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


from pysigmap.data import Data

# url = "".join(
#     [
#         "https://raw.githubusercontent.com/eamontoyaa/data4testing/",
#         "main/pysigmap/testData.csv",
#     ]
# )
# data = Data(
#     pd.read_csv(url),
#     sigmaV=75,
# )
# fig = data.plot()


url = "illustrative_examples/safet_mutapcija/Test_data_T11_simple.csv"

data = Data(
    pd.read_csv(url), sigmaV=75, reloading=False, secondUnloading=False
)
fig = data.plot()


# url = "illustrative_examples/safet_mutapcija/Test_data_T11.csv"

# data = Data(
#     pd.read_csv(url),
#     sigmaV=75,
#     use_this_UR_stage=1,
# )
# fig = data.plot()


# Default parameters
data.compressionIdx(range2fitCc=None)
data.recompressionIdx(opt=1)
fig = data.plot()


data.compressionIdx(range2fitCc=(150, np.inf))
data.recompressionIdx(opt=2)
fig = data.plot()


# data.compressionIdx(range2fitCc=(300, 5000))
# data.recompressionIdx(opt=3)
# fig = data.plot()


# data.idxCc = 0.25
# fig = data.plot()


# -- Casagrande method
from pysigmap.casagrande import Casagrande

method1 = Casagrande(data)
data.compressionIdx(range2fitCc=None)
fig1 = method1.getSigmaP()


# fig1 = method1.getSigmaP(range2fitFOP=[5, 5000], loglog=True)

# fig1 = method1.getSigmaP(range2fitFOP=[5, 5000], loglog=False)


# data.compressionIdx(range2fitCc=(700, np.inf))
# fig1 = method1.getSigmaP()


# data.compressionIdx(range2fitCc=(500, np.inf))
# fig1 = method1.getSigmaP(range2fitFOP=[20, np.inf], loglog=True)


# fig1 = method1.getSigmaP(mcp=200)
# fig1 = method1.getSigmaP()


# -- Becker et al. method
from pysigmap.energy import BeckerEtAl

method2 = BeckerEtAl(data)
data.compressionIdx(range2fitCc=None)
fig2 = method2.getSigmaP(zoom=1.2)


# fig2 = method2.getSigmaP(range2fitCR=[500, np.inf], zoom=3.5)


# data.compressionIdx(range2fitCc=(500, 4000))
# method2 = BeckerEtAl(data)
# fig2 = method2.getSigmaP(zoom=3.5)


# # -- Morin  method
# from pysigmap.energy import BeckerEtAl

# method3 = BeckerEtAl(data)
# fig3 = method3.getSigmaP(
#     range2fitCR=[500, np.inf], zoom=3.5, morinFormulation=True
# )


# -- Wang & Frost method
from pysigmap.energy import WangAndFrost

method4 = WangAndFrost(data)
data.compressionIdx(range2fitCc=None)
fig4 = method4.getSigmaP()


# fig4 = method4.getSigmaP(range2fitCR=[500, np.inf])


# data.compressionIdx(range2fitCc=(1000, 4000))
# fig4 = method4.getSigmaP()


# -- Bilogaritmic methods
from pysigmap.bilog import Bilog

method5 = Bilog(data)
data.compressionIdx(range2fitCc=None)
fig5_1 = method5.getSigmaP(opt=1)
# fig5_2 = method5.getSigmaP(opt=2)
# fig5_3 = method5.getSigmaP(opt=3)

# fig5_1 = method5.getSigmaP(
#     range2fitRR=[5, 100], range2fitCR=[1000, np.inf], opt=1
# )
# fig5_2 = method5.getSigmaP(
#     range2fitRR=[5, 30], range2fitCR=[1000, np.inf], opt=2
# )
# fig5_3 = method5.getSigmaP(
#     range2fitRR=[5, 30], range2fitCR=[1000, np.inf], opt=3
# )

# data.compressionIdx(range2fitCc=(1000, 4000))
# fig5_1 = method5.getSigmaP(opt=1)
# fig5_2 = method5.getSigmaP(opt=2)
# fig5_3 = method5.getSigmaP(opt=3)


# -- Pacheco-Silva's method
from pysigmap.pachecosilva import PachecoSilva

method6 = PachecoSilva(data)
data.compressionIdx(range2fitCc=None)
fig6 = method6.getSigmaP()

# data.compressionIdx(range2fitCc=(1000, np.inf))
# fig6 = method6.getSigmaP()

# -- Boone method
from pysigmap.boone import Boone

method7 = Boone(data)
data.compressionIdx(range2fitCc=None)
fig7 = method7.getSigmaP()

# data.compressionIdx(range2fitCc=(1000, np.inf))
# fig7 = method7.getSigmaP()

plt.show()
