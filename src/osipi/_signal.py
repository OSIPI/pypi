import numpy as np


def signal_linear(R1: np.float64, k: np.float64) -> np.float64:
    """Linear model for relationship between R1 and magnitude signal

    Args:
        R1 (np.float64): longitudinal relaxation rate in units of /s. [OSIPI code Q.EL1.001]
        k (np.float64): proportionality constant in a.u. S [OSIPI code Q.GE1.009]

    Returns:
        np.float64: magnitude signal in a.u. [OSIPI code Q.MS1.001]

    References:
        - Lexicon url: https://osipi.github.io/OSIPI_CAPLEX/perfusionModels/#LinModel_SM2
        - Lexicon code: M.SM2.001
        - OSIPI name: Linear model
        - Adapted from equation given in the Lexicon
    """
    # calculate signal
    return k * R1  # S


def signal_SPGR(R1: np.float64, S0: np.float64, TR: np.float64, a: np.float64) -> np.float64:
    """Steady-state signal for SPGR sequence.

    Args:
        R1 (np.float64): longitudinal relaxation rate in units of /s. [OSIPI code Q.EL1.001]
        S0 (np.float64): fully T1-relaxed signal in a.u. [OSIPI code Q.MS1.010]
        TR (np.float64): repetition time in units of s. [OSIPI code Q.MS1.006]
        a (np.float64): prescribed flip angle in units of deg. [OSIPI code Q.MS1.007]

    Returns:
        np.float64: magnitude signal in a.u. [OSIPI code Q.MS1.001]

    References:
        - Lexicon url: https://osipi.github.io/OSIPI_CAPLEX/perfusionModels/#SPGR%20model
        - Lexicon code: M.SM2.002
        - OSIPI name: Spoiled gradient recalled echo model
        - Adapted from equation given in the Lexicon and contribution from MJT_UoEdinburgh_UK
    """

    # calculate signal
    a_rad = a * np.pi / 180
    exp_TR_R1 = np.exp(-TR * R1)
    return S0 * (((1.0 - exp_TR_R1) * np.sin(a_rad)) / (1.0 - exp_TR_R1 * np.cos(a_rad)))  # S
