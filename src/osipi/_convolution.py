import numpy as np


def exp_conv(T: float, t: np.ndarray, a: np.ndarray) -> np.ndarray:
    """Exponential convolution operation of (1/T)exp(-t/T) with a

    Args:
        T (float): exponent in time units
        t (np.ndarray): array of time points
        a (np.ndarray): array to be convolved with time exponential

    Returns:
        np.ndarray: convolved array

    """

    if T == 0:
        return a

    n = len(t)
    f = np.zeros((n,))

    x = (t[1:n - 1] - t[0:n - 2]) / T
    da = (a[1:n - 1] - a[0:n - 2]) / x

    E = np.exp(-x)
    E0 = 1 - E
    E1 = x - E0

    add = a[0:n - 2] * E0 + da * E1

    for i in range(0, n - 2):
        f[i + 1] = E[i] * f[i] + add[i]

    f[n - 1] = f[n - 2]
    return f

