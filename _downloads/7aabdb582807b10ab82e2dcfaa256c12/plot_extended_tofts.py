"""
====================
The Extended Tofts model
====================

Simulating tissue concentrations from extended Tofts model with different settings.
"""

# %%
# Import necessary packages
import numpy as np
import matplotlib.pyplot as plt
import osipi

# %%
# Generate Parker AIF with default settings.

# Define time points in units of seconds - in this case we use a time resolution of 1 sec and a total duration of 6 minutes.
t = np.arange(0, 6*60, 1)

# Create an AIF with default settings
ca = osipi.aif_parker(t)

# %%
# Plot the tissue concentrations for an extracellular volume fraction of 0.2 and 3 different plasma volumes of 0.05, 0.2 and 0.6
Ktrans = 0.2 # in units of 1/min
ve = 0.2  # volume fraction between 0 and 1
vp = [0.05, 0.2, 0.6] # volume fraction between 0 and 1
ct = osipi.extended_tofts(t, ca, Ktrans, ve, vp[0])
plt.plot(t, ct, 'b-', label=f'vp = {vp[0]}')
ct = osipi.extended_tofts(t, ca, Ktrans, ve, vp[1])
plt.plot(t, ct, 'g-', label=f'vp = {vp[1]}')
ct = osipi.extended_tofts(t, ca, Ktrans, ve, vp[2])
plt.plot(t, ct, 'm-', label=f'vp = {vp[2]}')
plt.xlabel('Time (sec)')
plt.ylabel('Tissue concentration (mM)')
plt.legend()
plt.show()

# %%
# Comparing different discretization methods for an extracellular volume fraction of 0.2, Ktrans of 0.2 /min and vp of 0.05
ct = osipi.extended_tofts(t, ca, Ktrans, ve, vp[0]) # Defaults to Convolution
plt.plot(t, ct, 'b-', label='Convolution')
ct = osipi.extended_tofts(t, ca, Ktrans, ve, vp[0], discretization_method='exp')
plt.plot(t, ct, 'g-', label='Exponential Convolution')
plt.title(f'Ktrans = {Ktrans} /min')
plt.xlabel('Time (sec)')
plt.ylabel('Tissue concentration (mM)')
plt.legend()
plt.show()

# Choose the last image as a thumbnail for the gallery
# sphinx_gallery_thumbnail_number = -1
