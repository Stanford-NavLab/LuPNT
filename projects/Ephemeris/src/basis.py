import numpy as np


def get_fourier_basis(Phi, order):
    """
    Get the fourier basis
    """
    if isinstance(Phi, np.ndarray):
        lent = Phi.size
    else:
        lent = 1

    if order == 0:
        return None
    else:
        T = np.zeros((lent, order * 2))
        for i in range(order):
            T[:, 2 * i] = np.cos((i + 1) * 2 * Phi)
            T[:, 2 * i + 1] = np.sin((i + 1) * 2 * Phi)
        return T


def get_fourier_basis_dt(Phi, Phidot, order):
    """
    Get the fourier basis
    """
    if isinstance(Phi, np.ndarray):
        lent = Phi.size
    else:
        lent = 1

    if order == 0:
        return None
    else:
        T = np.zeros((lent, order * 2))
        for i in range(order):
            T[:, 2 * i] = -(i + 1) * 2 * Phidot * np.sin((i + 1) * 2 * Phi)
            T[:, 2 * i + 1] = (i + 1) * 2 * Phidot * np.cos((i + 1) * 2 * Phi)
        return T


def get_poly_basis(t_k, order, t_fit, poly_type):
    # pre-computation for chebyshev
    lent = t_k.size
    if poly_type == "monomial":
        Tp = np.zeros((lent, order + 1))
        Tp[:, 0] = 1
        for i in range(1, order + 1):
            Tp[:, i] = 1.0 / i * t_k**i
        return Tp

    elif poly_type == "chebyshev":
        z_k = 2 * t_k / t_fit  # for chebyshev  [-1, 1]
        # Chebyshev coefficients
        Tc = np.zeros((lent, order + 1))
        Tc[:, 0] = 1
        if order == 0:
            return Tc
        Tc[:, 1] = z_k
        for i in range(2, order + 1):
            Tc[:, i] = 2 * z_k * Tc[:, i - 1] - Tc[:, i - 2]
        return Tc
    elif poly_type == "legendre":
        z_k = 2 * t_k / t_fit
        # Legendre coefficients
        Pl = np.zeros((lent, order + 1))
        Pl[:, 0] = 1
        if order == 0:
            return Pl
        Pl[:, 1] = z_k
        for i in range(2, order + 1):
            Pl[:, i] = (2 * i - 1) / i * z_k * Pl[:, i - 1] - (i - 1) / i * Pl[:, i - 2]
        return Pl
    else:
        print("Invalid polynomial type: ", poly_type)
        return None


def get_poly_basis_dt(t_k, order, t_fit, poly_type):
    # pre-computation for chebyshev
    lent = t_k.size
    if poly_type == "monomial":
        Tpdot = np.zeros((lent, order + 1))
        Tpdot[:, 0] = 0
        Tpdot[:, 1] = 1
        for i in range(2, order + 1):
            Tpdot[:, i] = t_k ** (i - 1)
        return Tpdot

    elif poly_type == "chebyshev":
        z_k = 2 * t_k / t_fit
        Tc = get_poly_basis(t_k, order, t_fit, "chebyshev")
        # Chebyshev coefficients
        Tcdot = np.zeros((lent, order + 1))
        Tcdot[:, 0] = 0
        Tcdot[:, 1] = 2 / t_fit
        for i in range(2, order + 1):
            Tcdot[:, i] = (
                2 * z_k * Tcdot[:, i - 1] + 2 * Tc[:, i - 1] * (2 / t_fit) - Tcdot[:, i - 2]
            )

        return Tcdot

    elif poly_type == "legendre":
        z_k = 2 * t_k / t_fit
        Pl = get_poly_basis(t_k, order, t_fit, "legendre")
        # Legendre coefficients
        Pldot = np.zeros((lent, order + 1))
        Pldot[:, 0] = 0
        Pldot[:, 1] = 2 / t_fit
        for i in range(2, order + 1):
            Pldot[:, i] = (
                (2 * i - 1) / i * z_k * Pldot[:, i - 1]
                + (2 * i - 1) / i * Pl[:, i - 1] * (2 / t_fit)
                - (i - 1) / i * Pldot[:, i - 2]
            )

        return Pldot
