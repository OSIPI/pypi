import warnings

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import convolve

from ._convolution import exp_conv


def tofts(
    t: np.ndarray,
    ca: np.ndarray,
    Ktrans: float,
    ve: float,
    Ta: float = 30.0,
    discretization_method: str = "conv",
) -> np.ndarray:
    """Tofts model as defined by Tofts and Kermode (1991)

    Args:
        t (np.ndarray): array of time points in units of sec. [OSIPI code Q.GE1.004]
        ca (np.ndarray):
            Arterial concentrations in mM for each time point in t. [OSIPI code Q.IC1.001]
        Ktrans (float):
            Volume transfer constant in units of 1/min. [OSIPI code Q.PH1.008]
        ve (float):
            Relative volume fraction of the extracellular
            extravascular compartment (e). [OSIPI code Q.PH1.001.[e]]
        Ta (float, optional):
            Arterial delay time,
            i.e., difference in onset time between tissue curve and AIF in units of sec. Defaults to 30 seconds. [OSIPI code Q.PH1.007]
        discretization_method (str, optional): Defines the discretization method. Options include

            – 'conv': Numerical convolution (default) [OSIPI code G.DI1.001]

            – 'exp': Exponential convolution [OSIPI code G.DI1.006]


    Returns:
        np.ndarray: Tissue concentrations in mM for each time point in t.

    See Also:
        `extended_tofts`

    References:
        - Lexicon url:
            https://osipi.github.io/OSIPI_CAPLEX/perfusionModels/#indicator-kinetic-models
        - Lexicon code: M.IC1.004
        - OSIPI name: Tofts Model
        - Adapted from contributions by: LEK_UoEdinburgh_UK, ST_USyd_AUS, MJT_UoEdinburgh_UK

    Example:

        Create an array of time points covering 6 min in steps of 1 sec,
        calculate the Parker AIF at these time points, calculate tissue concentrations
        using the Tofts model and plot the results.

        Import packages:

        >>> import matplotlib.pyplot as plt
        >>> import osipi
        >>> import numpy

        Calculate AIF:

        >>> t = np.arange(0, 6 * 60, 1)
        >>> ca = osipi.aif_parker(t)

        Calculate tissue concentrations and plot:

        >>> Ktrans = 0.6  # in units of 1/min
        >>> ve = 0.2  # takes values from 0 to 1
        >>> ct = osipi.tofts(t, ca, Ktrans, ve)
        >>> plt.plot(t, ca, "r", t, ct, "b")

    """
    if not np.allclose(np.diff(t), np.diff(t)[0]):
        warnings.warn(
            ("Non-uniform time spacing detected. Time array may be" " resampled."),
            stacklevel=2,
        )

    if Ktrans <= 0 or ve <= 0:
        ct = 0 * ca

    else:
        # Convert units
        Ktrans = Ktrans / 60  # from 1/min to 1/sec

        if discretization_method == "exp":  # Use exponential convolution
            # Shift the AIF by the arterial delay time (if not zero)
            if Ta != 0:
                f = interp1d(
                    t,
                    ca,
                    kind="linear",
                    bounds_error=False,
                    fill_value=0,
                )
                ca = (t > Ta) * f(t - Ta)

            Tc = ve / Ktrans
            ct = ve * exp_conv(Tc, t, ca)

        else:  # Use convolution by default
            # Calculate the impulse response function
            kep = Ktrans / ve
            imp = Ktrans * np.exp(-1 * kep * t)

            # Shift the AIF by the arterial delay time (if not zero)
            if Ta != 0:
                f = interp1d(
                    t,
                    ca,
                    kind="linear",
                    bounds_error=False,
                    fill_value=0,
                )
                ca = (t > Ta) * f(t - Ta)

            # Check if time data grid is uniformly spaced
            if np.allclose(np.diff(t), np.diff(t)[0]):
                # Convolve impulse response with AIF
                convolution = np.convolve(ca, imp)

                # Discard unwanted points and make sure time spacing
                # is correct
                ct = convolution[0 : len(t)] * t[1]
            else:
                # Resample at the smallest spacing
                dt = np.min(np.diff(t))
                t_resampled = np.linspace(t[0], t[-1], int((t[-1] - t[0]) / dt))
                ca_func = interp1d(
                    t,
                    ca,
                    kind="quadratic",
                    bounds_error=False,
                    fill_value=0,
                )
                imp_func = interp1d(
                    t,
                    imp,
                    kind="quadratic",
                    bounds_error=False,
                    fill_value=0,
                )
                ca_resampled = ca_func(t_resampled)
                imp_resampled = imp_func(t_resampled)
                # Convolve impulse response with AIF
                convolution = np.convolve(ca_resampled, imp_resampled)

                # Discard unwanted points and make sure time spacing
                # is correct
                ct_resampled = convolution[0 : len(t_resampled)] * t_resampled[1]

                # Restore time grid spacing
                ct_func = interp1d(
                    t_resampled,
                    ct_resampled,
                    kind="quadratic",
                    bounds_error=False,
                    fill_value=0,
                )
                ct = ct_func(t)

    return ct


def extended_tofts(
    t: np.ndarray,
    ca: np.ndarray,
    Ktrans: float,
    ve: float,
    vp: float,
    Ta: float = 30.0,
    discretization_method: str = "conv",
) -> np.ndarray:
    """Extended tofts model as defined by Tofts (1997)

    Args:
        t (np.ndarray):
            array of time points in units of sec. [OSIPI code Q.GE1.004]
        ca (np.ndarray):
            Arterial concentrations in mM for each time point in t. [OSIPI code Q.IC1.001]
        Ktrans (float):
            Volume transfer constant in units of 1/min. [OSIPI code Q.PH1.008]
        ve (float):
            Relative volume fraction of the extracellular
            extravascular compartment (e). [OSIPI code Q.PH1.001.[e]]
        vp (float):
            Relative volyme fraction of the plasma compartment (p). [OSIPI code Q.PH1.001.[p]]
        Ta (float, optional):
            Arterial delay time, i.e., difference in onset time
            between tissue curve and AIF in units of sec.
            Defaults to 30 seconds. [OSIPI code Q.PH1.007]
        discretization_method (str, optional):
            Defines the discretization method. Options include

            – 'conv': Numerical convolution (default) [OSIPI code G.DI1.001]

            – 'exp': Exponential convolution [OSIPI code G.DI1.006]


    Returns:
        np.ndarray: Tissue concentrations in mM for each time point in t.

    See Also:
        `tofts`

    References:
        - Lexicon url: https://osipi.github.io/OSIPI_CAPLEX/perfusionModels/#indicator-kinetic-models
        - Lexicon code: M.IC1.005
        - OSIPI name: Extended Tofts Model
        - Adapted from contributions by: LEK_UoEdinburgh_UK, ST_USyd_AUS, MJT_UoEdinburgh_UK

    Example:

        Create an array of time points covering 6 min in steps of 1 sec,
        calculate the Parker AIF at these time points, calculate tissue concentrations
        using the Extended Tofts model and plot the results.

        Import packages:

        >>> import matplotlib.pyplot as plt
        >>> import osipi

        Calculate AIF

        >>> t = np.arange(0, 6 * 60, 0.1)
        >>> ca = osipi.aif_parker(t)

        Calculate tissue concentrations and plot

        >>> Ktrans = 0.6  # in units of 1/min
        >>> ve = 0.2  # takes values from 0 to 1
        >>> vp = 0.3  # takes values from 0 to 1
        >>> ct = osipi.extended_tofts(t, ca, Ktrans, ve, vp)
        >>> plt.plot(t, ca, "r", t, ct, "b")

    """

    if not np.allclose(np.diff(t), np.diff(t)[0]):
        warnings.warn(
            ("Non-uniform time spacing detected. Time array may be" " resampled."),
            stacklevel=2,
        )

    if Ktrans <= 0 or ve <= 0:
        ct = vp * ca

    else:
        # Convert units
        Ktrans = Ktrans / 60  # from 1/min to 1/sec

        if discretization_method == "exp":  # Use exponential convolution
            # Shift the AIF by the arterial delay time (if not zero)
            if Ta != 0:
                f = interp1d(
                    t,
                    ca,
                    kind="linear",
                    bounds_error=False,
                    fill_value=0,
                )
                ca = (t > Ta) * f(t - Ta)

            Tc = ve / Ktrans
            # expconv calculates convolution of ca and
            # (1/Tc)exp(-t/Tc), add vp*ca term for extended model
            ct = (vp * ca) + ve * exp_conv(Tc, t, ca)

        else:  # Use convolution by default
            # Calculate the impulse response function
            kep = Ktrans / ve
            imp = Ktrans * np.exp(-1 * kep * t)

            # Shift the AIF by the arterial delay time (if not zero)
            if Ta != 0:
                f = interp1d(
                    t,
                    ca,
                    kind="linear",
                    bounds_error=False,
                    fill_value=0,
                )
                ca = (t > Ta) * f(t - Ta)

            # Check if time data grid is uniformly spaced
            if np.allclose(np.diff(t), np.diff(t)[0]):
                # Convolve impulse response with AIF
                convolution = np.convolve(ca, imp)

                # Discard unwanted points, make sure time spacing is
                # correct and add vp*ca term for extended model
                ct = convolution[0 : len(t)] * t[1] + (vp * ca)
            else:
                # Resample at the smallest spacing
                dt = np.min(np.diff(t))
                t_resampled = np.linspace(t[0], t[-1], int((t[-1] - t[0]) / dt))
                ca_func = interp1d(
                    t,
                    ca,
                    kind="quadratic",
                    bounds_error=False,
                    fill_value=0,
                )
                imp_func = interp1d(
                    t,
                    imp,
                    kind="quadratic",
                    bounds_error=False,
                    fill_value=0,
                )
                ca_resampled = ca_func(t_resampled)
                imp_resampled = imp_func(t_resampled)
                # Convolve impulse response with AIF
                convolution = np.convolve(ca_resampled, imp_resampled)

                # Discard unwanted points, make sure time spacing is
                # correct and add vp*ca term for extended model
                ct_resampled = convolution[0 : len(t_resampled)] * t_resampled[1] + (
                    vp * ca_resampled
                )

                # Restore time grid spacing
                ct_func = interp1d(
                    t_resampled,
                    ct_resampled,
                    kind="quadratic",
                    bounds_error=False,
                    fill_value=0,
                )
                ct = ct_func(t)

    return ct


def two_compartment_exchange_model(
    t: np.ndarray,
    ca: np.ndarray,
    Fp: float,
    Ps: float,
    Ve: float,
    Vp: float,
    Ta: float = 30.0,
) -> np.ndarray:
    """Two compartment exchange model

    Tracer flows from the AIF to the blood plasma compartment; two-way leakage
    between the plasma and extracellular compartments(EES) is permitted.

    Args:
        t: 1D array of times(s).
        ca: Arterial concentrations in mM for each time point in t.
        Fp: Blood plasma flow rate into a unit tissue volume in ml/min.
        Ps: Permeability surface area product in ml/min.
        Ve: Extracellular volume fraction.
        Vp: Plasma volume fraction.
        Ta: Arterial delay time, i.e., difference in onset time between tissue curve and AIF in units of sec.

    Returns:
        Ct: Tissue concentrations in mM

    References:
        - Lexicon url: https://osipi.github.io/OSIPI_CAPLEX/perfusionModels/#indicator-kinetic-models
        - Lexicon code: M.IC1.009
        - OSIPI name: Two Compartment Exchange Model
        - Adapted from contributions by: MJT_UoEdinburgh_UK

    """
    if not np.allclose(np.diff(t), np.diff(t)[0]):
        warnings.warn(
            ("Non-uniform time spacing detected. Time array may be" " resampled."),
            stacklevel=2,
        )

    # Convert units
    fp_per_s = Fp / (60.0 * 100.0)
    ps_per_s = Ps / 60.0

    # Calculate the impulse response function
    v = Ve + Vp

    # Mean transit time
    T = v / fp_per_s
    tc = Vp / fp_per_s
    te = Ve / ps_per_s

    sig_p = ((T + te) + np.sqrt((T + te) ** 2 - 4 * tc * te)) / (2 * tc * te)
    sig_n = ((T + tc) - np.sqrt((T + tc) ** 2 - 4 * tc * te)) / (2 * tc * te)

    # Calculate the impulse response function for the plasma compartment and EES

    irf_cp = (
        Vp
        * sig_p
        * sig_n
        * ((1 - te * sig_n) * np.exp(-t * sig_n) + (te * sig_p - 1.0) * np.exp(-t * sig_p))
        / (sig_p - sig_n)
    )

    irf_ce = Ve * sig_p * sig_n * (np.exp(-t * sig_n) - np.exp(-t * sig_p)) / (sig_p - sig_n)

    irf_cp[[0]] /= 2
    irf_ce[[0]] /= 2

    dt = np.min(np.diff(t))

    # get concentration in plasma and EES

    Cp = dt * convolve(ca, irf_cp, mode="full", method="auto")[: len(t)]
    Ce = dt * convolve(ca, irf_ce, mode="full", method="auto")[: len(t)]

    # get tissue concentration
    Ct = Cp + Ce
    return Ct
