"""
==============
Signal to concentration via electromagnetic property (SPGR, FXL, analytical linear relaxivity)
==============

Demonstrating the SPGR model for relationship between
    signal S and total tissue indicator concentration, assuming the FXL.
"""

# %%
# Import necessary packages
import matplotlib.pyplot as plt
import numpy as np
import osipi

# %%
# Convert a series of S values to the corresponding indicator concentrations.
# Example data adapted from OSIPI repository ("vox_1")
S = np.array(
    [
        7,
        9,
        6,
        10,
        9,
        6,
        9,
        10,
        9,
        9,
        9,
        12,
        8,
        10,
        12,
        15,
        53,
        70,
        71,
        70,
        63,
        58,
        54,
        50,
        53,
        48,
        52,
        49,
        42,
    ],
    dtype=np.float64,
)
S_baseline = S[0]  # use first point for baseline signal
R10 = np.float64(1 / 1.4)
TR = np.float64(0.002)
a = np.float64(13)
r1 = np.float64(4.5)

C = osipi.S_to_C_via_R1_SPGR(S, S_baseline, R10, TR, a, r1)
print(f"Concentration (mM): {C}")

# Plot S and C
fig, ax = plt.subplots(2, 1)
ax[0].plot(S, "b-")
ax[0].set_ylabel("S (a.u.)")
ax[1].plot(C, "b-")
ax[1].set_ylabel("C (mM)")
ax[1].set_xlabel("time point")
plt.show()
