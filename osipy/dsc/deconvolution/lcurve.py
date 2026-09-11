"""L-curve selection of the Tikhonov regularization parameter.

DSC deconvolution recovers R(t) from C(t) = CBF * AIF(t) ** R(t). The AIF
convolution matrix is severely ill-conditioned, so the minimum-residual
solution amplifies noise instead of recovering R(t). Tikhonov regularization
damps that by replacing the raw inverse 1/s_i with s_i / (s_i^2 + lam^2).

The L-curve criterion picks lambda from the data instead of a hard-coded
constant like in (TSVD). Plotting solution norm against residual norm over a range of lambda
traces an L-shaped curve; its corner sits where the fit stops explaining
signal and starts explaining noise, and is found as the point of maximum
curvature.

All functions take the array module as ``xp`` and work on numpy or cupy.

References
----------
.. [1] Hansen PC, O'Leary DP. The use of the L-curve in the regularization
   of discrete ill-posed problems. SIAM J Sci Comput 1993;14(6):1487-1503.
   doi:10.1137/0914086
.. [2] Calamante F, Gadian DG, Connelly A. Quantification of bolus-tracking
   MRI: improved characterization of the tissue residue function using
   Tikhonov regularization. MRM 2003;50(6):1237-1247. doi:10.1002/mrm.10643
.. [3] OSIPI DCE-DSC-MRI CodeCollection, SR_TBG_BNIPhoenix_USA
   AIFDeconvolution, https://github.com/OSIPI/DCE-DSC-MRI_CodeCollection
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray

#: Number of points on the logarithmic lambda grid.
DEFAULT_N_LAMBDA = 400

#: Grid bounds as a fraction of the largest singular value. Scaling by
#: s_max keeps the search invariant to concentration units and time step.
DEFAULT_LAMBDA_MIN = 1e-9
DEFAULT_LAMBDA_MAX = 10.0

#: Floor for denominators that vanish at both ends of the grid.
_EPS = 1e-30


def select_lambda(
    S: NDArray[np.floating[Any]],
    UtC: NDArray[np.floating[Any]],
    xp: Any,
    n_lambda: int = DEFAULT_N_LAMBDA,
    lambda_min: float = DEFAULT_LAMBDA_MIN,
    lambda_max: float = DEFAULT_LAMBDA_MAX,
) -> NDArray[np.floating[Any]]:
    """Pick a per-voxel regularization parameter at the L-curve corner.

    With b = U.T @ C the data projected onto the left singular vectors, the
    Tikhonov solution for a given lambda has squared residual norm, squared
    solution norm, and derivative of the solution norm:

        rho  = sum_i [ lam^2 / (s_i^2 + lam^2) ]^2 * b_i^2
        eta  = sum_i [ s_i   / (s_i^2 + lam^2) ]^2 * b_i^2
        eta' = sum_i -4 * lam * s_i^2 / (s_i^2 + lam^2)^3 * b_i^2

    The curvature of (log rho, log eta) is then::

        g = 2 * (eta * rho / eta')
              * (lam^2 * eta' * rho + 2 * lam * eta * rho
                 + lam^4 * eta * eta')
              / (lam^4 * eta^2 + rho^2)^1.5

    Following Hansen this is the *negative* of the curvature, so the corner
    is a minimum of g, not a maximum.

    Parameters
    ----------
    S : NDArray
        Singular values of the AIF convolution matrix, shape ``(n_sv,)``.
    UtC : NDArray
        ``U.T @ C``, shape ``(n_sv, n_voxels)``.
    xp : module
        Array module (``numpy`` or ``cupy``).
    n_lambda : int
        Number of grid points.
    lambda_min, lambda_max : float
        Grid bounds relative to the largest singular value.

    Returns
    -------
    NDArray
        Chosen lambda per voxel, shape ``(n_voxels,)``.

    """
    s_max = float(xp.max(S))
    if s_max <= 0.0:
        return xp.zeros(UtC.shape[1], dtype=UtC.dtype)

    # Log-spaced grid, scaled by s_max so the bounds are relative.
    k = xp.arange(n_lambda, dtype=UtC.dtype) / max(n_lambda - 1, 1)
    lam = s_max * lambda_min * (lambda_max / lambda_min) ** k  # (n_lambda,)

    # Weights are shared across voxels; only b^2 varies.
    lam_c = lam[:, None]
    denom = S[None, :] ** 2 + lam_c**2  # (n_lambda, n_sv)
    w_rho = (lam_c**2 / denom) ** 2
    w_eta = (S[None, :] / denom) ** 2
    w_deta = -4.0 * lam_c * S[None, :] ** 2 / denom**3

    b2 = UtC**2  # (n_sv, n_voxels)
    rho = w_rho @ b2  # (n_lambda, n_voxels)
    eta = w_eta @ b2
    eta_prime = w_deta @ b2

    # eta' and the denominator both vanish at the ends of the grid, where
    # there is no corner to find. Floor them and zero g there.
    num = (
        lam_c**2 * eta_prime * rho
        + 2.0 * lam_c * eta * rho
        + lam_c**4 * eta * eta_prime
    )
    den = (lam_c**4 * eta**2 + rho**2) ** 1.5
    valid = (xp.abs(eta_prime) > _EPS) & (den > _EPS)

    g = xp.where(
        valid,
        2.0
        * (eta * rho / xp.where(valid, eta_prime, -1.0))
        * num
        / xp.where(valid, den, 1.0),
        0.0,
    )

    return lam[_last_local_min(g, xp)]


def _last_local_min(g: NDArray[np.floating[Any]], xp: Any) -> NDArray[np.integer[Any]]:
    """Index of the highest-lambda local minimum of g, per voxel.

    Equivalent to the reference implementation's scan down from the top of
    the grid, stopping at the first strict local minimum. Voxels with no
    interior local minimum fall back to the global minimum.

    Parameters
    ----------
    g : NDArray
        Negative curvature, shape ``(n_lambda, n_voxels)``.
    xp : module
        Array module (``numpy`` or ``cupy``).

    Returns
    -------
    NDArray
        Grid index per voxel, shape ``(n_voxels,)``.
    """
    interior = g[1:-1]
    is_min = (interior < g[:-2]) & (interior < g[2:])

    # argmax on the reversed axis gives the last True, so the scan runs
    # from large lambda downwards.
    n_interior = is_min.shape[0]
    last = n_interior - 1 - xp.argmax(is_min[::-1], axis=0) + 1
    return xp.where(xp.any(is_min, axis=0), last, xp.argmin(g, axis=0))


def tikhonov_filter_factors(
    S: NDArray[np.floating[Any]],
    lambdas: NDArray[np.floating[Any]],
    xp: Any,
) -> NDArray[np.floating[Any]]:
    """Tikhonov filter factors s_i / (s_i^2 + lam^2).

    Replaces the raw pseudo-inverse 1/s_i. Truncated SVD applies a hard
    0/1 cutoff; this rolls off smoothly instead, damping modes near the
    cutoff rather than dropping them.

    Parameters
    ----------
    S : NDArray
        Singular values, shape ``(n_sv,)``.
    lambdas : NDArray
        Lambda per voxel, shape ``(n_voxels,)``.
    xp : module
        Array module (``numpy`` or ``cupy``).

    Returns
    -------
    NDArray
        Filter factors, shape ``(n_sv, n_voxels)``.
    """
    denom = S[:, None] ** 2 + lambdas[None, :] ** 2
    return xp.where(denom > 0.0, S[:, None] / xp.where(denom > 0.0, denom, 1.0), 0.0)


def solve_tikhonov_lcurve(
    U: NDArray[np.floating[Any]],
    S: NDArray[np.floating[Any]],
    Vh: NDArray[np.floating[Any]],
    concentration: NDArray[np.floating[Any]],
    xp: Any,
    n_lambda: int = DEFAULT_N_LAMBDA,
    lambda_min: float = DEFAULT_LAMBDA_MIN,
    lambda_max: float = DEFAULT_LAMBDA_MAX,
) -> tuple[NDArray[np.floating[Any]], NDArray[np.floating[Any]]]:
    """Deconvolve with Tikhonov filtering, choosing lambda by the L-curve.

    Convenience entry point shared by both deconvolution paths.

    Parameters
    ----------
    U, S, Vh : NDArray
        SVD of the AIF convolution matrix.
    concentration : NDArray
        Tissue concentration curves, shape ``(n_timepoints, n_voxels)``.
    xp : module
        Array module (``numpy`` or ``cupy``).
    n_lambda : int
        Number of grid points.
    lambda_min, lambda_max : float
        Grid bounds relative to the largest singular value.

    Returns
    -------
    residue : NDArray
        Residue functions, shape ``(n_rows_of_Vh, n_voxels)``.
    lambdas : NDArray
        Lambda chosen per voxel, shape ``(n_voxels,)``.
    """
    UtC = U.T @ concentration
    lambdas = select_lambda(
        S, UtC, xp, n_lambda=n_lambda, lambda_min=lambda_min, lambda_max=lambda_max
    )
    residue = Vh.T @ (tikhonov_filter_factors(S, lambdas, xp) * UtC)
    return residue, lambdas
