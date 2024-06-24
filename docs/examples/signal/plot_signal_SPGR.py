"""
==============
Linear model for relationship between R1 and magnitude signal
==============

Demonstrating the linear model for relationship between R1 and magnitude signal, s = k.R1
"""

# %%
# Import necessary packages
import matplotlib.pyplot as plt
import numpy as np
import osipi

# %%
# Convert a series of R1 values to the corresponding signal intensities using the SPGR model.

R1 = np.linspace(0.1, 10, 50)  # R1 in units of /s.
S0 = np.float64(100)  # fully T1-relaxed signal in a.u.
TR = np.float64(5e-3)  # repetition time in units of s.
a = np.float64(15)  # prescribed flip angle in units of deg.
S = osipi.signal_SPGR(R1, S0, TR, a)  # signal in a.u.
print(f"Signal: {S}")

# Plot S vs. R1
plt.plot(R1, S, "r-")
plt.xlabel("R1 (/s)")
plt.ylabel("S (a.u.)")
plt.show()
