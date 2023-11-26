import pylupnt as pnt
import numpy as np
import sys
import pytest

sys.path.append("pylupntutil")
import pntautodiff as ad


def Amat():
    A = np.zeros((4, 3))
    for i in range(4):
        for j in range(3):
            A[i, j] = 1.0 / (i + j * 2 + 1)

    return A


def myfunc(x):
    A = Amat()
    y = ad.realvec(A @ x, "X")
    return y


def myfunc2(x):
    y = ad.divvecscaler(ad.expvec(x), ad.sumvec(x))
    return y


def test_ad_jacobian_1():
    x_val = [0.5, 0.2, 0.3]
    x = ad.realvec(x_val, "X")
    A = Amat()
    J = pnt.jacobian(myfunc, x)
    np.testing.assert_array_almost_equal(J, A)


def test_ad_jacobian_2():
    x_val = [0.5, 0.2, 0.3]
    x = ad.realvec(x_val, "X")
    J = pnt.jacobian(myfunc2, x)
    Jnum = ad.numerical_jacobian(myfunc2, x)
    np.testing.assert_array_almost_equal(J, Jnum)


if __name__ == "__main__":
    test_ad_jacobian_1()
    test_ad_jacobian_2()
    print("All tests passed")
