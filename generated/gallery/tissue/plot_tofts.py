"""
====================
The Tofts model
====================

Simulating tissue concentrations from Tofts model with different settings.
"""

import matplotlib.pyplot as plt

# %%
# Import necessary packages
import numpy as np
import osipi

# %%
# Generate Parker AIF with default settings.

# Define time points in units of seconds - in this case we use a time
# resolution of 1 sec and a total duration of 6 minutes.
t = np.arange(0, 6 * 60, 1)

# Create an AIF with default settings
ca = osipi.aif_parker(t)

# %%
# Plot the tissue concentrations for an extracellular volume fraction
# of 0.2 and 3 different transfer rate constants of 0.05, 0.2 and 0.6
# /min
Ktrans = [0.05, 0.2, 0.6]  # in units of 1/min
ve = 0.2  # volume fraction between 0 and 1
ct = osipi.tofts(t, ca, Ktrans=Ktrans[0], ve=ve)
plt.plot(t, ct, "b-", label=f"Ktrans = {Ktrans[0]} /min")
ct = osipi.tofts(t, ca, Ktrans[1], ve)
plt.plot(t, ct, "g-", label=f"Ktrans = {Ktrans[1]} /min")
ct = osipi.tofts(t, ca, Ktrans[2], ve)
plt.plot(t, ct, "m-", label=f"Ktrans = {Ktrans[2]} /min")
plt.xlabel("Time (sec)")
plt.ylabel("Tissue concentration (mM)")
plt.legend()
plt.show()

# %%
# Comparing different discretization methods for an extracellular
# volume fraction of 0.2 and Ktrans of 0.2 /min
ct = osipi.tofts(t, ca, Ktrans=Ktrans[1], ve=ve)  # Defaults to Convolution
plt.plot(t, ct, "b-", label="Convolution")
ct = osipi.tofts(t, ca, Ktrans=Ktrans[1], ve=ve, discretization_method="exp")
plt.plot(t, ct, "g-", label="Exponential Convolution")
plt.title(f"Ktrans = {Ktrans[1]} /min")
plt.xlabel("Time (sec)")
plt.ylabel("Tissue concentration (mM)")
plt.legend()
plt.show()

# Choose the last image as a thumbnail for the gallery
# sphinx_gallery_thumbnail_number = -1
