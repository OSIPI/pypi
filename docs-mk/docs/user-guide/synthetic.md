## Read MR data from a folder of dicom files
``` python
import osipi


dicom_folder = 'path/to/dicom/folder'

signal, slices, dicom_ref = osipi.read_dicom(dicom_folder)

#Read a DICOM series from a folder path.
#Returns the signal, slices and dicom reference.
# - The signal is a 4D numpy array with dimensions (x, y, z, t)
# - The slices is a list of slices
# - The dicom reference is a sample dicom file from the folder

```

## Enhance the signal
Enhance the signal by removing baseline signal which is the first time point of the signal.
Before the contrast agent is arrived, the signal is assumed to be constant. This constant signal is removed from the signal to enhance the signal.
``` python
import osipi
from osipi import enhance_signal

E, S0, S = enhance_signal(signal, data_shape, 5)

# - E is the enhanced signal after removing the baseline signal
# - S0 is the average baseline signal here in this example of the first 5 time points
# - S is the original raw signal
```

## Get the AIF
Get the Arterial Input Function (AIF) from the signal.
Using a mask, the AIF is extracted from the signal.

``` python
import osipi
from osipi import get_aif_from_ROI
# first you may have to create a mask for the AIF manually or using a saved mask and apply it to the signal
aif_mask, rio_voxels = osipi.rio(signal, slice, saved=True)
# this returns the mask and the number of voxels in the mask
aif = get_aif_from_ROI(signal, aif_mask)
# this returns the AIF from the signal by averaging the signal over the voxels in the mask
```

## Get the perfusion and tissue parameters

Get the perfusion and tissue parameters from the signal and the AIF.
You may here choose different models to fit the signal to get the parameters.

``` python
import osipi
from osipi import extended_tofts_model

# Fit the signal to the extended Tofts model
ktrans, ve, vp = extended_tofts_model(signal, aif, data_shape)

# - ktrans is the volume transfer constant
# - ve is the extravascular extracellular volume fraction
# - vp is the plasma volume fraction
# visit CAPLEX for more information on the model and the parameters
```
