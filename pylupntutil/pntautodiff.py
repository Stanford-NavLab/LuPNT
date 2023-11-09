import numpy as np
import pylupnt as pnt


# initializations
def realvec_with_size(xval, n=None):
    if isinstance(xval, np.ndarray):
        size = xval.size
    else:
        size = len(xval)

    if n == "X":
        x = pnt.VectorXreal(size)
    else:
        x = getattr(pnt, f"Vector{size}real")(size)

    for i in range(size):
        x[i] = xval[i]
    return x

# in default, vectors are dynamic
def realvec(xval):
    return realvec_with_size(xval, n="X")

def realvecn(xval, n):
    return realvec_with_size(xval, n=n)

# vector functions
def funcvec(x, f):
    n = len(x)
    z = pnt.VectorXreal(n)
    for i in range(n):
        z[i] = f(x[i])
    return z

def normvec(x):
    return pnt.sqrt(dotvec(x, x))


def absvec(x):
    return funcvec(x, lambda t: pnt.abs(t))


def expvec(x):
    return funcvec(x, lambda t: pnt.exp(t))


def logvec(x):
    return funcvec(x, lambda t: pnt.log(t))


def log10vec(x):
    return funcvec(x, lambda t: pnt.log10(t))


def cbrtvec(x):
    return funcvec(x, lambda t: pnt.cbrt(t))


def sinvec(x):
    return funcvec(x, lambda t: pnt.sin(t))


def cosvec(x):
    return funcvec(x, lambda t: pnt.cos(t))


def tanvec(x):
    return funcvec(x, lambda t: pnt.tan(t))


def arcsinvec(x):
    return funcvec(x, lambda t: pnt.arcsin(t))


def arccosvec(x):
    return funcvec(x, lambda t: pnt.arccos(t))


def arctanvec(x):
    return funcvec(x, lambda t: pnt.arctan(t))


def sinhvec(x):
    return funcvec(x, lambda t: pnt.sinh(t))


def coshvec(x):
    return funcvec(x, lambda t: pnt.cosh(t))


def tanhvec(x):
    return funcvec(x, lambda t: pnt.tanh(t))


def arcsinhvec(x):
    return funcvec(x, lambda t: pnt.arcsinh(t))


def arccoshvec(x):
    return funcvec(x, lambda t: pnt.arccosh(t))


def sumvec(x):
    n = len(x)
    s = pnt.real(0.0)
    for i in range(n):
        s += x[i]
    return s


def arctan2vec(x, y):
    n = len(x)
    z = pnt.VectorXreal(n)
    for i in range(n):
        z[i] = pnt.arctan2(x[i], y[i])
    return z


def maxvec(x, y):
    n = len(x)
    z = pnt.VectorXreal(n)
    for i in range(n):
        z[i] = pnt.max(x[i], y[i])
    return z


def minvec(x, y):
    n = len(x)
    z = pnt.VectorXreal(n)
    for i in range(n):
        z[i] = pnt.min(x[i], y[i])
    return z


def powvec(x, y):
    n = len(x)
    z = pnt.VectorXreal(n)
    for i in range(n):
        z[i] = pnt.pow(x[i], y[i])
    return z


def dotvec(x, y):
    n = len(x)
    z = pnt.real(0.0)
    for i in range(n):
        z += x[i] * y[i]
    return z


def multvec(x, y):
    n = len(x)
    z = pnt.VectorXreal(n)
    for i in range(n):
        z[i] = x[i] * y[i]
    return z


def divvec(x, y):
    n = len(x)
    z = pnt.VectorXreal(n)
    for i in range(n):
        z[i] = x[i] / y[i]
    return z


def divvecscaler(x, y):
    n = len(x)
    z = pnt.VectorXreal(n)
    for i in range(n):
        z[i] = x[i] / y
    return z


def numerical_gradient(f, x, eps=1e-8):
    n = len(x)
    grad_num = np.zeros(n)
    x_val = x.asarray()

    for i in range(n):
        x_p = pnt.VectorXreal(n)
        x_m = pnt.VectorXreal(n)
        for j in range(n):
            if i == j:
                x_p[j] = x_val[j] + eps
                x_m[j] = x_val[j] - eps
            else:
                x_p[j] = x_val[j]
                x_m[j] = x_val[j]
        F_eps_p = f(x_p)
        F_eps_m = f(x_m)
        grad_num[i] = (F_eps_p - F_eps_m) / (2.0 * eps)

    return grad_num


# matrix functions
def numerical_jacobian(f, x, eps=1e-8):
    n = len(x)
    m = len(f(x))
    jac_num = np.zeros((m, n))
    x_val = x.asarray()

    for i in range(n):
        x_p = pnt.VectorXreal(n)
        x_m = pnt.VectorXreal(n)
        for j in range(n):
            if i == j:
                eps_ = max(abs(x_val[j] * eps), eps)
                x_p[j] = x_val[j] + eps_
                x_m[j] = x_val[j] - eps_
            else:
                x_p[j] = x_val[j]
                x_m[j] = x_val[j]
        F_eps_p = f(x_p)
        F_eps_m = f(x_m)
        jac_num[:, i] = (F_eps_p.asarray() - F_eps_m.asarray()) / (2.0 * eps_)

    return jac_num
