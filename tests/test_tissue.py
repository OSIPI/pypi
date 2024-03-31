import numpy as np
import osipi
import math
import matplotlib.pyplot as plt
import pytest


@pytest.mark.parametrize("method", ["tofts", "extended_tofts"])
def test_tissue_concentration_model(method):
    if method == "tofts":
        function = osipi.tofts
    elif method == "extended_tofts":
        function = osipi.extended_tofts

    # 1. Basic operation of the function - test that the peak tissue concentration is less than the peak AIF
    t = np.linspace(0, 6 * 60, 360)
    ca = osipi.aif_parker(t)
    ct = function(t, ca, Ktrans=0.6, ve=0.2, vp=0.3 if method == "extended_tofts" else None)
    assert np.round(np.max(ct)) < np.round(np.max(ca))

    # 2. Basic operation of the function - test with non-uniform spacing of time array
    t = np.geomspace(1, 6*60+1, num=360)-1
    ca = osipi.aif_parker(t)
    ct = function(t, ca, Ktrans=0.6, ve=0.2, vp=0.3 if method == "extended_tofts" else None)
    assert np.round(np.max(ct)) < np.round(np.max(ca))

    # 3. The offset option - test that the tissue concentration is shifted from the AIF by the specified offset time
    t = np.arange(0, 6 * 60, 1)
    ca = osipi.aif_parker(t)
    ct = function(t, ca, Ktrans=0.6, ve=0.2, vp=0.3 if method == "extended_tofts" else None, Ta=60.0)
    assert (np.min(np.where(ct > 0.0)) - np.min(np.where(ca > 0.0)) - 1) * 1 == 60.0

    # 4. Test that the discretization options give almost the same result - time step must be very small
    t = np.arange(0, 6 * 60, 0.01)
    ca = osipi.aif_parker(t)
    ct_conv = function(t, ca, Ktrans=0.6, ve=0.2, vp=0.3 if method == "extended_tofts" else None)
    ct_exp = function(t, ca, Ktrans=0.6, ve=0.2, vp=0.3 if method == "extended_tofts" else None, discretization_method='exp')
    assert np.allclose(ct_conv, ct_exp, rtol=1e-4, atol=1e-3)

    # 5. Test that the ratio of the area under the ct and ca curves is approximately the extracellular volume plus the plasma volume
    t = np.arange(0, 6 * 60, 1)
    ca = osipi.aif_parker(t)
    ct_conv = function(t, ca, Ktrans=0.6, ve=0.2, vp=0.3 if method == "extended_tofts" else None)
    ct_exp = function(t, ca, Ktrans=0.6, ve=0.2, vp=0.3 if method == "extended_tofts" else None, discretization_method='exp')
    assert math.isclose(np.trapz(ct_conv, t)/np.trapz(ca, t), 0.2 + 0.3 if method == "extended_tofts" else 0, abs_tol=1e-1)
    assert math.isclose(np.trapz(ct_exp, t)/np.trapz(ca, t), 0.2 + 0.3 if method == "extended_tofts" else 0, abs_tol=1e-1)

    # 6. Test specific use cases
    t = np.arange(0, 6 * 60, 1)
    ca = osipi.aif_parker(t)

    if function == "tofts":
        # Test case 1: Ktrans=0, ve=0.2
        ct_conv = osipi.tofts(t, ca, Ktrans=0, ve=0.2)
        assert np.count_nonzero(ct_conv) == 0

        ct_exp = osipi.tofts(t, ca, Ktrans=0, ve=0.2, discretization_method='exp')
        assert np.count_nonzero(ct_exp) == 0

        # Test case 2: Ktrans=0.6, ve=0
        ct_conv = osipi.tofts(t, ca, Ktrans=0.6, ve=0)
        assert np.count_nonzero(ct_conv) == 0

        ct_exp = osipi.tofts(t, ca, Ktrans=0.6, ve=0, discretization_method='exp')
        assert np.count_nonzero(ct_exp) == 0

    elif function == "extended_tofts":
        # Test case 1: Ktrans=0, ve=0.2, vp=0.3
        ct_conv = osipi.extended_tofts(t, ca, Ktrans=0, ve=0.2, vp=0.3)
        assert np.allclose(ct_conv, ca * 0.3, rtol=1e-4, atol=1e-3)

        ct_exp = osipi.extended_tofts(t, ca, Ktrans=0, ve=0.2, vp=0.3, discretization_method='exp')
        assert np.allclose(ct_conv, ca * 0.3, rtol=1e-4, atol=1e-3)

        # Test case 2: Ktrans=0.6, ve=0, vp=0.3
        ct_conv = osipi.extended_tofts(t, ca, Ktrans=0.6, ve=0, vp=0.3)
        assert np.allclose(ct_conv, ca * 0.3, rtol=1e-4, atol=1e-3)

        ct_exp = osipi.extended_tofts(t, ca, Ktrans=0.6, ve=0, vp=0.3, discretization_method='exp')
        assert np.allclose(ct_conv, ca * 0.3, rtol=1e-4, atol=1e-3)


if __name__ == "__main__":

    test_tissue_concentration_model()

    print('All tissue concentration model tests passed!!')

