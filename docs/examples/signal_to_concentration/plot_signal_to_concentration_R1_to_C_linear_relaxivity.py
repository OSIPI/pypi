"""
==============
Electromagnetic property inverse model: longitudinal relaxation rate, linear with relaxivity
==============

Demonstrating the inverse linear relaxivity model
    for converting R_1 to tissue indicator concentration C.
"""

# %%
# Import necessary packages
import matplotlib.pyplot as plt
import numpy as np
import osipi

# %%
# Convert a series of R1 values to the corresponding indicator concentrations.
R1 = np.array([1, 2, 3, 4, 5, 6], dtype=np.float64)
R10 = np.float64(1)
r1 = np.float64(5)

C = osipi.R1_to_C_linear_relaxivity(R1, R10, r1)
print(f"Concentration (mM): {C}")

# Plot C vs. R1
plt.plot(R1, C, "r-")
plt.xlabel("R1 (/s)")
plt.ylabel("C (mM)")
plt.show()
