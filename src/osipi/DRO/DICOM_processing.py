import os

import numpy as np
import pydicom

from osipi.DRO.filters_and_noise import median_filter


def read_dicom_slices_as_4d_signal(folder_path: str) -> (np.ndarray, dict, np.ndarray):
    """Read DICOM files from a folder and return the 4D signal data.

    Args:
        folder_path (str): Path to the folder containing DICOM files.

    Returns:
        signal: The 4D signal data (x, y, z, t).
        slices: Dictionary containing the slices.
        example_data: The first slice.


    """

    slices = {}
    for root, _, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".dcm"):
                dicom_file = os.path.join(root, file)
                slice = pydicom.read_file(dicom_file)
                if slice.SliceLocation not in slices:
                    slices[slice.SliceLocation] = []
                slices[slice.SliceLocation].append((slice.AcquisitionTime, slice))

    # Sort each list of slices by the first element (AcquisitionTime)
    for slice_location in slices:
        slices[slice_location].sort(key=lambda x: x[0])

    spatial_shape = slices[slice_location][0][1].pixel_array.shape

    data_shape = (spatial_shape[0], spatial_shape[1], len(slices), len(slices[slice_location]))

    signal = np.zeros(data_shape)

    for z, slice_location in enumerate(sorted(slices.keys())):  # Sort by slice location
        for t, (_, slice) in enumerate(
            sorted(slices[slice_location], key=lambda x: x[0])
        ):  # Sort by acquisition time
            signal[:, :, z, t] = slice.pixel_array

    return signal, slices, slices[slice_location][0][1]


def SignalEnhancementExtract(
    S: np.ndarray, datashape: list, baselinepoints: int = 5
) -> (np.ndarray, np.ndarray, np.ndarray):
    """Signal enhancement extraction from 4D signal data.

    Args:
        S: is the 4D signal data (x, y, z, t).
        datashape: is the shape of the data (x, y, z, t).
        baselinepoints: is the number of time points before contrast injection.

    Returns:
        E: The signal enhancement.
        S0: The baseline signal (S0).
        S: The 4D signal data (x, y, z, t).

    """
    S0 = np.average(S[:, :, :, 0:baselinepoints], axis=3)  # Take baseline signal
    E = np.zeros_like(S)

    # Calcualte siganl enhancement
    for i in range(0, datashape[-1]):
        E[:, :, :, i] = S[:, :, :, i] - S0
        E[:, :, :, i] = median_filter(E[:, :, :, i])  # Median filter size (3,3)

    return E, S0, S


def calculate_baseline(signal, baseline):
    """
    Calculate the baseline signal (S0) from pre-contrast time points.

    Parameters:
    signal (numpy.ndarray): The 4D signal data (x, y, z, t).
    baseline (int): Number of time points before contrast injection.

    Returns:
    numpy.ndarray: The baseline signal (S0).
    """
    S0 = np.average(signal[:, :, :, :baseline], axis=3, keepdims=True)
    return S0


def signal_to_R1(signal, S0, TR):
    """
    Convert signal to R1 values using the Ernst equation.

    Parameters:
    signal (numpy.ndarray): The 4D signal data (x, y, z, t).
    S0 (numpy.ndarray): The baseline signal (S0).
    TR (float): Repetition time.

    Returns:
    numpy.ndarray: The R1 values.
    """
    epsilon = 1e-8  # Small constant to avoid division by zero and log of zero
    R1 = -1 / TR * np.log((signal + epsilon) / (S0 + epsilon))
    return R1


def calc_concentration(R1, R10, r1):
    """
    Calculate the concentration of the contrast agent in tissue.

    Parameters:
    R1 (numpy.ndarray): The R1 values.
    R10 (numpy.ndarray): The pre-contrast R1 values.
    r1 (float): Relaxivity of the contrast agent.

    Returns:
    numpy.ndarray: The concentration of the contrast agent in the tissue.
    """
    Ctiss = (R1 - R10) / r1
    return Ctiss


def signal_enhancement(signal, S0, R10, r1):
    """
    Calculate the signal enhancement.

    Parameters:
    signal (numpy.ndarray): The 4D signal data (x, y, z, t).
    other parameters same as previous function

    Returns:
    numpy.ndarray: The signal enhancement.
    """
    E = (R10 / r1) * (signal - S0) / S0
    return E
