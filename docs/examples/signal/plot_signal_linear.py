"""
==============
Linear model for relationship between R_1 and magnitude signal S
==============

Demonstrating the linear model for relationship between R_1 and magnitude signal, S = k.R_1
"""

# %%
# Import necessary packages
import matplotlib.pyplot as plt
import numpy as np
import osipi

# %%
# Convert a series of R1 values to the corresponding signal intensities.

R1 = np.array([0.0, 1.5, 3.0, 4.0, 10.0])  # R_1 in units of /s
k = np.float64(150.0)  # constant of proportionality in units of a.u. s
S = osipi.signal_linear(R1, k)  # signal in a.u.
print(f"Signal: {S}")

# Plot S vs. R1
plt.plot(R1, S, "ro-")
plt.xlabel("R1 (/s)")
plt.ylabel("S (a.u.)")
plt.show()
