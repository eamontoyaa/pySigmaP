=====
Usage
=====

Loading input data from an external file
----------------------------------------

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    fig = data.plot()
    plt.show()


Computing the compression and recompression indices
---------------------------------------------------

* Default parameters: :math:`C_\mathrm{c}` (maximum slope) and :math:`C_\mathrm{r}` (opt=1)

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    data.compressionIdx(range2fitCc=None)
    data.recompressionIdx(opt=1)
    fig = data.plot()
    plt.show()

* :math:`C_\mathrm{c}` (two last points) and :math:`C_\mathrm{r}` (opt=2)

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    data.compressionIdx(range2fitCc=(3000, np.inf))
    data.recompressionIdx(opt=2)
    fig = data.plot()
    plt.show()

* :math:`C_\mathrm{c}` (four last points) and :math:`C_\mathrm{r}` (opt=3)

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    data.compressionIdx(range2fitCc=(700, np.inf))
    data.recompressionIdx(opt=2)
    fig = data.plot()
    plt.show()

Casagrande method
-----------------

* Default parameters: cubic spline function

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data
    from pysigmap.casagrande import Casagrande

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    method = Casagrande(data)
    fig = method.getSigmaP(mcp=None, range2fitFOP=None, loglog=True)
    plt.show()

* Fourth order polynomial (FOP)

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data
    from pysigmap.casagrande import Casagrande

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    method = Casagrande(data)
    fig = method.getSigmaP(range2fitFOP=[20, 5000], loglog=True)
    plt.show()

* Maximum curvature point (MCP) manually introduced

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data
    from pysigmap.casagrande import Casagrande

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    method = Casagrande(data)
    fig = method.getSigmaP(mcp=200)
    plt.show()

Pacheco Silva method
--------------------

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data
    from pysigmap.pachecosilva import PachecoSilva

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    method = PachecoSilva(data)
    fig = method.getSigmaP()
    plt.show()

Boone method
------------

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data
    from pysigmap.boone import Boone

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    method = Boone(data)
    fig = method.getSigmaP()
    plt.show()


Bilogarithmic methods
---------------------

* Butterfield method

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data
    from pysigmap.bilog import Bilog

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    method = Bilog(data)
    fig = method.getSigmaP(range2fitRR=None, range2fitCR=None, opt=1)
    plt.show()

* Oikawa method

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data
    from pysigmap.bilog import Bilog

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    method = Bilog(data)
    fig = method.getSigmaP(range2fitRR=None, range2fitCR=[1000, 5000], opt=2)
    plt.show()

* Onitsuka et al. method

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data
    from pysigmap.bilog import Bilog

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    method = Bilog(data)
    fig = method.getSigmaP(range2fitRR=[0, 30], range2fitCR=[1000, 9000], opt=3)
    plt.show()


Strain energy methods
---------------------

* Becker et al. method

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data
    from pysigmap.energy import BeckerEtAl

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    method = BeckerEtAl(data)
    fig = method.getSigmaP(range2fitRR=None, range2fitCR=None, morinFormulation=False, zoom=5.5)
    plt.show()

* Morin method

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data
    from pysigmap.energy import BeckerEtAl

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    method = BeckerEtAl(data)
    fig = method.getSigmaP(range2fitRR=[0, 100], range2fitCR=[700, 9000], morinFormulation=True, zoom=5.5)
    plt.show()

* Wang and Frost method

.. plot::
    :include-source:
    :align: center

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pysigmap.data import Data
    from pysigmap.energy import WangAndFrost

    url = "https://raw.githubusercontent.com/eamontoyaa/data4testing/main/pysigmap/testData.csv"
    data = Data(pd.read_csv(url), sigmaV=75)
    method = WangAndFrost(data)
    fig = method.getSigmaP(range2fitCR=None)
    plt.show()
