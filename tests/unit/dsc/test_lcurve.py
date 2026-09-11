"""Unit tests for L-curve selection of the Tikhonov regularization parameter."""

import numpy as np
import pytest

from osipy.dsc.deconvolution.lcurve import (
    _last_local_min,
    select_lambda,
    solve_tikhonov_lcurve,
    tikhonov_filter_factors,
)


def build_problem(
    n_timepoints: int = 60,
    dt: float = 1.0,
    cbf: float = 50.0,
    mtt: float = 4.0,
    noise_frac: float = 0.01,
    seed: int = 0,
):
    """Build a deconvolution problem with a known residue function.

    Returns the SVD of the convolution matrix plus the noisy data and the
    ground truth, so tests can check recovery rather than only shapes.
    """
    rng = np.random.default_rng(seed)
    time = np.arange(n_timepoints) * dt

    shifted = time - 5.0
    aif = np.where(shifted > 0, (shifted**3) * np.exp(-shifted / 1.5), 0.0)
    aif = aif / aif.max() * 6.0

    residue = cbf * np.exp(-time / mtt)

    idx = np.arange(n_timepoints)[:, None] - np.arange(n_timepoints)[None, :]
    A = np.where(idx >= 0, aif[idx], 0.0) * dt

    clean = A @ residue
    noisy = clean + rng.normal(0.0, noise_frac * clean.max(), n_timepoints)

    U, S, Vh = np.linalg.svd(A)
    return {
        "A": A,
        "U": U,
        "S": S,
        "Vh": Vh,
        "cbf_true": cbf,
        "noisy": noisy[:, None],
        "noise_norm": float(np.linalg.norm(noisy - clean)),
    }


class TestSelectLambda:
    """Behaviour of the corner selection."""

    def test_shape_and_positivity(self):
        p = build_problem()
        lam = select_lambda(p["S"], p["U"].T @ p["noisy"], np)
        assert lam.shape == (1,)
        assert lam[0] > 0

    def test_lands_in_a_useful_range(self):
        """Guards the corner rule itself.

        Picking the wrong extremum of g puts lambda at either end of the
        grid: a global maximum gives lambda/s_max ~ 0.8 (CBF collapses to
        ~11) and a global minimum gives ~1e-6 (no regularization at all).
        """
        p = build_problem()
        lam = select_lambda(p["S"], p["U"].T @ p["noisy"], np)
        assert 1e-5 < lam[0] / p["S"].max() < 0.3

    def test_more_noise_selects_more_regularization(self):
        """The corner has to move with the noise level, not sit still."""
        lams = []
        for noise_frac in (0.002, 0.005, 0.01, 0.02, 0.05):
            p = build_problem(noise_frac=noise_frac, seed=1)
            lams.append(select_lambda(p["S"], p["U"].T @ p["noisy"], np)[0])

        assert np.all(np.diff(lams) > 0), lams

    def test_invariant_to_concentration_scaling(self):
        """Lambda lives in the singular-value space of A, not the data's.

        Scaling the data by alpha scales rho, eta and eta_prime by
        alpha**2, so g's prefactor goes as alpha**2, its numerator as
        alpha**4 and its denominator as alpha**6 -- alpha cancels.
        """
        p = build_problem()
        UtC = p["U"].T @ p["noisy"]

        lam = select_lambda(p["S"], UtC, np)
        lam_scaled = select_lambda(p["S"], 1000.0 * UtC, np)
        assert lam_scaled[0] == pytest.approx(lam[0], rel=1e-9)

    def test_tracks_scaling_of_the_convolution_matrix(self):
        """Scaling A (units, dt) must scale lambda with it."""
        p = build_problem()
        UtC = p["U"].T @ p["noisy"]

        lam = select_lambda(p["S"], UtC, np)
        lam_scaled = select_lambda(100.0 * p["S"], 100.0 * UtC, np)
        assert lam_scaled[0] == pytest.approx(100.0 * lam[0], rel=1e-6)

    def test_voxels_are_independent(self):
        """Column j of a batch must equal the same voxel run on its own."""
        cols = [
            build_problem(noise_frac=f, seed=s)["noisy"][:, 0]
            for f, s in ((0.002, 1), (0.02, 2), (0.05, 3))
        ]
        p = build_problem()
        batch = np.stack(cols, axis=1)

        lam_batch = select_lambda(p["S"], p["U"].T @ batch, np)
        assert lam_batch.shape == (3,)

        for j in range(3):
            alone = select_lambda(p["S"], p["U"].T @ batch[:, j : j + 1], np)
            assert lam_batch[j] == pytest.approx(alone[0], rel=1e-12)

    def test_handles_all_zero_singular_values(self):
        S = np.zeros(5)
        UtC = np.ones((5, 2))
        lam = select_lambda(S, UtC, np)
        assert lam.shape == (2,)
        assert np.all(np.isfinite(lam))

    def test_handles_zero_data(self):
        p = build_problem()
        lam = select_lambda(p["S"], np.zeros((p["S"].size, 1)), np)
        assert np.all(np.isfinite(lam))


class TestLastLocalMin:
    """The corner rule: last local minimum, not the first and not global."""

    def test_picks_the_highest_index_local_minimum(self):
        g = np.array([5.0, 1.0, 5.0, 4.0, 2.0, 4.0, 9.0])[:, None]
        assert int(_last_local_min(g, np)[0]) == 4

    def test_ignores_a_deeper_minimum_at_lower_index(self):
        """A global argmin would return 1 here; the scan must return 4."""
        g = np.array([5.0, -99.0, 5.0, 4.0, 2.0, 4.0, 9.0])[:, None]
        assert int(np.argmin(g[:, 0])) == 1
        assert int(_last_local_min(g, np)[0]) == 4

    def test_falls_back_to_global_min_when_monotonic(self):
        g = np.arange(6.0)[::-1][:, None]
        assert int(_last_local_min(g, np)[0]) == 5

    def test_per_voxel(self):
        g = np.array(
            [
                [5.0, 9.0],
                [1.0, 8.0],
                [5.0, 1.0],
                [4.0, 7.0],
                [2.0, 6.0],
                [4.0, 5.0],
            ]
        )
        assert list(_last_local_min(g, np)) == [4, 2]


class TestTikhonovFilterFactors:
    """Filter factor properties."""

    def test_shape(self):
        filt = tikhonov_filter_factors(
            np.array([4.0, 2.0, 1.0]), np.array([0.1, 0.5]), np
        )
        assert filt.shape == (3, 2)

    def test_approaches_inverse_for_small_lambda(self):
        S = np.array([4.0, 2.0, 1.0])
        filt = tikhonov_filter_factors(S, np.array([1e-8]), np)
        assert np.allclose(filt[:, 0], 1.0 / S, rtol=1e-6)

    def test_damps_modes_below_lambda(self):
        """A mode well below lambda must be suppressed, not amplified."""
        S = np.array([10.0, 1.0, 1e-4])
        filt = tikhonov_filter_factors(S, np.array([0.5]), np)
        assert filt[2, 0] < 1.0
        assert filt[2, 0] < 1e-6 * (1.0 / S[2])

    def test_zero_singular_value_gives_zero_not_nan(self):
        filt = tikhonov_filter_factors(np.array([1.0, 0.0]), np.array([0.1]), np)
        assert np.all(np.isfinite(filt))
        assert filt[1, 0] == pytest.approx(0.0)

    def test_zero_lambda_and_zero_singular_value(self):
        filt = tikhonov_filter_factors(np.array([0.0]), np.array([0.0]), np)
        assert np.all(np.isfinite(filt))


class TestSolveTikhonovLcurve:
    """End-to-end deconvolution."""

    def test_shapes(self):
        p = build_problem()
        residue, lam = solve_tikhonov_lcurve(p["U"], p["S"], p["Vh"], p["noisy"], np)
        assert residue.shape == (60, 1)
        assert lam.shape == (1,)

    def test_matches_manual_composition(self):
        p = build_problem()
        residue, lam = solve_tikhonov_lcurve(p["U"], p["S"], p["Vh"], p["noisy"], np)

        UtC = p["U"].T @ p["noisy"]
        expected = p["Vh"].T @ (tikhonov_filter_factors(p["S"], lam, np) * UtC)
        assert np.allclose(residue, expected)

    def test_recovers_cbf_without_tuning(self):
        """The point of the method: no hand-set knob, CBF still comes back."""
        p = build_problem(cbf=50.0, mtt=4.0)
        residue, _ = solve_tikhonov_lcurve(p["U"], p["S"], p["Vh"], p["noisy"], np)
        assert residue[:, 0].max() == pytest.approx(50.0, rel=0.15)

    def test_recovers_cbf_across_noise_levels(self):
        for noise_frac in (0.002, 0.01, 0.05):
            p = build_problem(noise_frac=noise_frac, seed=1)
            residue, _ = solve_tikhonov_lcurve(p["U"], p["S"], p["Vh"], p["noisy"], np)
            assert residue[:, 0].max() == pytest.approx(50.0, rel=0.15), noise_frac

    def test_residual_lands_near_the_noise_floor(self):
        """Under-regularized fits drive the residual below the noise norm."""
        p = build_problem()
        residue, _ = solve_tikhonov_lcurve(p["U"], p["S"], p["Vh"], p["noisy"], np)

        resid = np.linalg.norm(p["A"] @ residue[:, 0] - p["noisy"][:, 0])
        assert 0.2 * p["noise_norm"] < resid < 2.0 * p["noise_norm"]

    def test_solution_does_not_blow_up(self):
        """The unregularized solve of this same problem reaches ~1e15."""
        p = build_problem()
        residue, _ = solve_tikhonov_lcurve(p["U"], p["S"], p["Vh"], p["noisy"], np)
        assert np.abs(residue).max() < 10.0 * p["cbf_true"]

    def test_residue_decays(self):
        """R(t) should start high and end near zero."""
        p = build_problem()
        residue, _ = solve_tikhonov_lcurve(p["U"], p["S"], p["Vh"], p["noisy"], np)
        r = residue[:, 0]
        assert r[0] > 0.5 * r.max()
        assert abs(r[-1]) < 0.1 * r.max()
