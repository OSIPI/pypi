import numpy as np
import osipi
import math

def test_tissue_tofts():
    Ktrans = 0.6
    ve = 0.2
    t = np.arange(0, 6 * 60, 0.01)
    ca = osipi.aif_parker(t)

    # 1. Basic operation of the function - test that the peak tissue concentration is less than the peak AIF
    ct = osipi.tofts(t, ca, Ktrans=Ktrans, ve=ve)
    assert np.round(np.max(ct)) < np.round(np.max(ca))

    # 2. The offset option - test that the tissue concentration is shifted from the AIF by the specified offset time
    ct = osipi.tofts(t, ca, Ktrans=Ktrans, ve=ve, t_offset=60.0)
    assert (np.min(np.where(ct>0.0)) - np.min(np.where(ca>0.0)) - 1)*0.01 == 60.0

    # 3. Test that the discretization options give almost the same result
    ct_conv = osipi.tofts(t, ca, Ktrans=0.6, ve=0.2)
    ct_exp = osipi.tofts(t, ca, Ktrans=0.6, ve=0.2, discretization_method='exp')
    assert np.allclose(ct_conv, ct_exp, rtol=1e-4, atol=1e-3)

    # 4. Test that the ratio of the area under the ct and ca curves is approximately the extracellular volume
    ct_conv = osipi.tofts(t, ca, Ktrans=Ktrans, ve=ve)
    ct_exp = osipi.tofts(t, ca, Ktrans=Ktrans, ve=ve, discretization_method='exp')
    assert math.isclose(np.trapz(ct_conv, t)/np.trapz(ca, t), ve, abs_tol=1e-1)
    assert math.isclose(np.trapz(ct_exp, t)/np.trapz(ca, t), ve, abs_tol=1e-1)

def test_tissue_extended_tofts():

    # Not implemented yet so need to raise an error
    Ktrans = 0.6
    ve = 0.2
    t = np.arange(0, 6 * 60, 0.01)
    ca = osipi.aif_parker(t)
    try:
        ct = osipi.extended_tofts(t, ca, Ktrans=Ktrans, ve=ve)
    except:
        assert True
    else:
        assert False

if __name__ == "__main__":

    test_tissue_tofts()
    test_tissue_extended_tofts()

    print('All tissue concentration model tests passed!!')

