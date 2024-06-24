
Signal to relaxation time
^^^^^^^^^^^^^^^^^^^^^^^^^

Convert a series of T1-weighted signal values to the corresponding R1 values.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import osipi

    S = np.array([7, 9, 6, 10, 9, 6, 9, 10, 9, 9, 9, 12, 8, 10, 12, 15, 53, 70, 71, 70, 63, 58, 54, 50, 53, 48, 52, 49, 42], dtype=np.float64)
    S_baseline = np.float64(S[0])  # use first point for baseline signal
    R10 = np.float64(1/1.4)
    TR = np.float64(.002)
    a = np.float64(13)
    R1 = osipi.S_to_R1_SPGR(S, S_baseline, R10, TR, a)
    plt.plot(R1)
    plt.ylabel('R1 (/s)')
    plt.show()

Relaxation time to concentration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Convert a series of R1 values to the corresponding indicator concentrations.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import osipi

    R1 = np.array([1, 2, 3, 4, 5, 6], dtype=np.float64)
    R10 = np.float64(1)
    r1 = np.float64(5)

    C = osipi.R1_to_C_linear_relaxivity(R1, R10, r1)
    plt.plot(C)
    plt.ylabel('C (mM)')
    plt.show()

Concentration to tissue parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Coming soon..

All in one go: signal to tissue parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Coming soon..
