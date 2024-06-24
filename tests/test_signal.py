import numpy as np
import osipi


def test_signal_linear():
    # 1. Simple use case to generate a single signal value
    R1 = 3.0  # R1 in units of /s
    k = np.float64(150.0)  # constant of proportionality in units of arb. unit s
    S_truth = 450.0  # expected signal in arb. unit
    S = osipi.signal_linear(R1, k)  # estimated signal
    np.testing.assert_allclose(S_truth, S, rtol=0, atol=1e-7)

    # 2. Simple use case to generate an array of signal values
    R1 = np.array([0.0, 1.5, 3.0, 4.0, 10.0], dtype=np.float64)  # R1 in units of /s
    k = np.float64(150.0)  # constant of proportionality in units of arb. unit s
    S_truth = np.array([0.0, 225.0, 450.0, 600.0, 1500.0])  # expected signal in arb. unit
    S = osipi.signal_linear(R1, k)  # estimated signal
    np.testing.assert_allclose(S_truth, S, rtol=0, atol=1e-7)


def test_signal_SPGR():
    # 1. Simple use case to generate a single signal value
    # test data created using https://github.com/mjt320/SEPAL
    R1 = 0.5  # R1 in units of /s
    S0 = np.float64(100)
    TR = np.float64(5e-3)
    a = np.float64(15)
    S_truth = 1.77119982
    S = osipi.signal_SPGR(R1, S0, TR, a)
    np.testing.assert_allclose(S_truth, S, rtol=0, atol=1e-7)

    # 2. Simple use case to generate an array of signal values
    # test data created using https://github.com/mjt320/SEPAL
    R1 = np.array([0.1, 0.2, 0.5, 1, 2, 10, 50], dtype=np.float64)  # R1 in units of /s
    S0 = np.float64(100)
    TR = np.float64(5e-3)
    a = np.float64(15)
    S_truth = np.array(
        [0.37438758, 0.73827771, 1.77119982, 3.31912401, 5.89510144, 15.5485315, 23.10948814],
        dtype=np.float64,
    )
    S = osipi.signal_SPGR(R1, S0, TR, a)
    np.testing.assert_allclose(S_truth, S, rtol=0, atol=1e-7)


if __name__ == "__main__":
    test_signal_linear()
    test_signal_SPGR()

    print("All signal functionality tests passed!!")
