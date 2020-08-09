"""test.py module."""


def main():
    """Script for testing the package."""
    import pandas as pd
    from pysimgap.data import Data

    # data = Data(pd.read_csv('testData1.csv'), sigmaV=100)
    data = Data(pd.read_csv('testDataCRB2.csv'), sigmaV=100, strainPercent=True,
                reloadStage=False, doubleUnload=False)

    # -- Casagrande method
    from pysimgap.casagrande import Casagrande
    method1 = Casagrande(data)
    fig1 = method1.getSigmaP(range2fitTOP=[25, 800], range2fitNCL=[500, 2000])

    # -- Becker et al. method
    from pysimgap.energy import BeckerEtAl
    method2 = BeckerEtAl(data)
    fig2 = method2.getSigmaP(range2fitNCL=[500, 2000], zoom=2.3)

    # -- Morin  method
    from pysimgap.energy import BeckerEtAl
    method3 = BeckerEtAl(data, morinApproach=True)
    fig3 = method3.getSigmaP(range2fitNCL=[500, 2000], zoom=2.3)

    # -- Wang & Frost method
    from pysimgap.energy import WangAndFrost
    method4 = WangAndFrost(data)
    fig4 = method4.getSigmaP(range2fitNCL=[500, 2000])

    # -- Bilogaritmic methods
    from pysimgap.bilog import Bilog
    method5 = Bilog(data)
    fig5_1 = method5.getSigmaP(range2fitNCL=[500, 2000], opt=1)
    fig5_2 = method5.getSigmaP(range2fitNCL=[500, 2000], opt=2)
    fig5_3 = method5.getSigmaP(range2fitNCL=[500, 2000], opt=3)

    # -- Pacheco-Silva's method
    from pysimgap.pachecosilva import PachecoSilva
    method6 = PachecoSilva(data)
    fig6 = method6.getSigmaP(range2fitTOP=[100, 400], range2fitNCL=[500, 2000])

    # -- Boone method
    from pysimgap.boone import Boone
    method7 = Boone(data)
    fig7 = method7.getSigmaP(range2fitTOP=[25, 800], range2fitNCL=[500, 2000])

if __name__ == '__main__':
    main()
