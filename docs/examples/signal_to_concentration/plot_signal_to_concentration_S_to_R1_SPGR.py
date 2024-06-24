"""
==============
Signal to electromagnetic properties (R1) (SPGR, FXL, analytical linear relaxivity)
==============

Demonstrating the SPGR model for relationship between signal S and R1, assuming the FXL.
"""

# %%
# Import necessary packages
import matplotlib.pyplot as plt
import numpy as np
import osipi

# %%
# Convert a series of S values to the corresponding R1 values.
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
S_baseline = np.float64(S[0])  # use first point for baseline signal
R10 = np.float64(1 / 1.4)
TR = np.float64(0.002)
a = np.float64(13)

R1 = osipi.S_to_R1_SPGR(S, S_baseline, R10, TR, a)
print(f"R_1 (/s): {R1}")

# Plot S and R1
fig, ax = plt.subplots(2, 1)
ax[0].plot(S, "b-")
ax[0].set_ylabel("S (a.u.)")
ax[1].plot(R1, "b-")
ax[1].set_ylabel("R1 (/s)")
ax[1].set_xlabel("time point")
plt.show()
