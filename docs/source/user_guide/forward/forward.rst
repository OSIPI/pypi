
Generating an AIF
^^^^^^^^^^^^^^^^^

Generate an AIF and plot it:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import osipi

    t = np.arange(0, 6*60, 1)
    ca = osipi.aif_parker(t)
    plt.plot(t, ca)
    plt.show()


Generating a tissue concentration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generate tissue concentration and plot it:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import osipi

    t = np.arange(0, 6*60, 1)
    ca = osipi.aif_parker(t)
    Ktrans = 0.6
    ve = 0.2
    ct = osipi.tofts(t, ca, Ktrans=Ktrans, ve=ve)
    plt.plot(t, ct)
    plt.show()


Generating an MRI signal
^^^^^^^^^^^^^^^^^^^^^^^^

Coming soon..


Adding measurement error
^^^^^^^^^^^^^^^^^^^^^^^^

Coming soon..


