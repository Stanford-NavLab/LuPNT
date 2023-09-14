import pylupnt as pnt
import numpy as np
import sys

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
    y = ad.realvec(A @ x)
    return y


def myfunc2(x):
    y = ad.divvecscaler(ad.expvec(x), ad.sumvec(x))
    return y


# main script -----------------------------------------------

x_val = [0.5, 0.2, 0.3]
x = ad.realvec(x_val)
print("x = ")
print(x)
print(" ")

print("F = Ax")
print("----------------------------")
# The A matrix
print("A = ")
print(Amat())
print(" ")

F = myfunc(x)
print("F = ")
print(F)
print(" ")

J = pnt.jacobian(myfunc, x)
print("Jacobian (autodiff)= ")
print(J)
print(" ")

Jnum = ad.numerical_jacobian(myfunc, x)
print("Jacobian (numerical) = ")
print(Jnum)

print(" ")
print("F = exp(x) / sum(x)")
print("----------------------------")
F = myfunc2(x)
print("F = ")
print(F)
print(" ")

J = pnt.jacobian(myfunc2, x)
print("Jacobian (autodiff)= ")
print(J)
print(" ")

Jnum = ad.numerical_jacobian(myfunc2, x)
print("Jacobian (numerical) = ")
print(Jnum)
