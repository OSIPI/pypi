import numpy as np
import osipi


def test_S_to_C_via_R1_SPGR():
    # 1. Simple use case
    # Test data and ground truth adapted from OSIPI repository tests ("vox_1")
    S = np.array(
        [
            7,
            9,
            6,
            10,
            9,
            6,
            9,
            10,
            9,
            9,
            9,
            12,
            8,
            10,
            12,
            15,
            53,
            70,
            71,
            70,
            63,
            58,
            54,
            50,
            53,
            48,
            52,
            49,
            42,
        ],
        dtype=np.float64,
    )
    S_baseline = S[0]  # use first point for baseline signal
    R10 = np.float64(1 / 1.4)
    TR = np.float64(0.002)
    a = np.float64(13)
    r1 = np.float64(4.5)

    C = osipi.S_to_C_via_R1_SPGR(S, S_baseline, R10, TR, a, r1)
    C_truth = np.array(
        [
            5.896517841898054e-15,
            0.04861114658132921,
            -0.02373686966428807,
            0.0735037381246513,
            0.04861114658132921,
            -0.02373686966428807,
            0.04861114658132921,
            0.0735037381246513,
            0.04861114658132921,
            0.04861114658132921,
            0.04861114658132921,
            0.12451099240893969,
            0.02411300135228399,
            0.0735037381246513,
            0.12451099240893969,
            0.2042308274298612,
            1.723816850668738,
            2.9860012006598615,
            3.0813949356951134,
            2.9860012006598615,
            2.3932805439961378,
            2.0365845677672656,
            1.7832505654332091,
            1.5539858164348437,
            1.723816850668738,
            1.4473321937838983,
            1.6658291725704821,
            1.5000331327922662,
            1.15513218792702,
        ],
        dtype=np.float64,
    )
    np.testing.assert_allclose(C_truth, C, rtol=0, atol=1e-7)


def test_S_to_R1_SPGR():
    # 1. Simple use case
    # Test data adapted from OSIPI repository tests ("vox_1")
    S = np.array(
        [
            7,
            9,
            6,
            10,
            9,
            6,
            9,
            10,
            9,
            9,
            9,
            12,
            8,
            10,
            12,
            15,
            53,
            70,
            71,
            70,
            63,
            58,
            54,
            50,
            53,
            48,
            52,
            49,
            42,
        ],
        dtype=np.float64,
    )
    S_baseline = S[0]  # use first point for baseline signal
    R10 = np.float64(1 / 1.4)
    TR = np.float64(0.002)
    a = np.float64(13)

    R1 = osipi.S_to_R1_SPGR(S, S_baseline, R10, TR, a)
    R1_truth = np.array(
        [
            0.71428571,
            0.93303587,
            0.6074698,
            1.04505254,
            0.93303587,
            0.6074698,
            0.93303587,
            1.04505254,
            0.93303587,
            0.93303587,
            0.93303587,
            1.27458518,
            0.82279422,
            1.04505254,
            1.27458518,
            1.63332444,
            8.47146154,
            14.15129112,
            14.58056292,
            14.15129112,
            11.48404816,
            9.87891627,
            8.73891326,
            7.70722189,
            8.47146154,
            7.22728059,
            8.21051699,
            7.46443481,
            5.91238056,
        ],
        dtype=np.float64,
    )
    np.testing.assert_allclose(R1_truth, R1, rtol=0, atol=1e-7)
    return


def test_R1_to_C_linear_relaxivity():
    # 1. Simple use case
    R1 = np.array([1, 2, 3, 4, 5, 6], dtype=np.float64)
    R10 = np.float64(1)
    r1 = np.float64(5)

    C = osipi.R1_to_C_linear_relaxivity(R1, R10, r1)
    C_truth = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], dtype=np.float64)
    np.testing.assert_allclose(C_truth, C, rtol=0, atol=1e-7)
    return


if __name__ == "__main__":
    test_S_to_C_via_R1_SPGR()
    test_S_to_R1_SPGR()
    test_R1_to_C_linear_relaxivity()

    print("All signal-to-concentration functionality tests passed!!")
