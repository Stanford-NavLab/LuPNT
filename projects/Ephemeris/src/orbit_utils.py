import numpy as np
import pylupnt as pnt


def wrapToPi(x):
    """
    Wrap the angle to the range [-pi, pi]

    x: float or array

    return: float or array
    """
    x = np.arctan2(np.sin(x), np.cos(x))
    return x


def wrapTo2Pi(x):
    """
    Wrap the angle to the range [0, 2*pi]
    """
    x = np.arctan2(np.sin(x), np.cos(x))  # wrap to [-pi, pi]
    x = x + 2 * np.pi * (x < 0)  # wrap to [0, 2*pi]
    return x


def rms(x):
    """
    Compute the root mean square of the input array.

    x: array

    return: float
    """
    return np.sqrt(np.mean(x**2))


def RotX(theta):
    return np.array(
        [[1, 0, 0], [0, np.cos(theta), np.sin(theta)], [0, -np.sin(theta), np.cos(theta)]]
    )


def RotY(theta):
    return np.array(
        [[np.cos(theta), 0, -np.sin(theta)], [0, 1, 0], [np.sin(theta), 0, np.cos(theta)]]
    )


def RotZ(theta):
    return np.array(
        [[np.cos(theta), np.sin(theta), 0], [-np.sin(theta), np.cos(theta), 0], [0, 0, 1]]
    )


def RotX_dot(theta, theta_dot):
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array(
        [[0, 0, 0], [0, -theta_dot * s, theta_dot * c], [0, -theta_dot * c, -theta_dot * s]]
    )


def RotY_dot(theta, theta_dot):
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array(
        [[-theta_dot * s, 0, -theta_dot * c], [0, 0, 0], [theta_dot * c, 0, -theta_dot * s]]
    )


def RotZ_dot(theta, theta_dot):
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array(
        [[-theta_dot * s, theta_dot * c, 0], [-theta_dot * c, -theta_dot * s, 0], [0, 0, 0]]
    )


def mean2ecc(M, e):
    """
    Convert mean anomaly to eccentric anomaly using Newton's method.

    Parameters
    ----------
    M : float
        Mean anomaly (radians).

    e : float
        Eccentricity.

    Returns
    -------
    E : float
        Eccentric anomaly (radians).
    """
    E = M + e * np.sin(M)  # Initial guess
    for _ in range(10):  # Iterate to refine the guess
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E


def ecc2true(E, e):
    """
    Convert eccentric anomaly to true anomaly.

    Parameters
    ----------
    E : float
        Eccentric anomaly (radians).

    e : float
        Eccentricity.

    Returns
    -------
    nu : float
        True anomaly (radians).
    """
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    return nu


def ecc2mean(E, e):
    """
    Convert eccentric anomaly to mean anomaly.

    Parameters
    ----------
    E : float
        Eccentric anomaly (radians).

    e : float
        Eccentricity.

    Returns
    -------
    M : float
        Mean anomaly (radians).
    """
    M = E - e * np.sin(E)
    return M


def mean_to_true(M, e):
    """
    Convert mean anomaly to true anomaly.

    Parameters
    ----------
    M : float
        Mean anomaly (radians).

    e : float
        Eccentricity.

    Returns
    -------
    nu : float
        True anomaly (radians).
    """
    E = mean2ecc(M, e)
    nu = ecc2true(E, e)
    return nu


def true2ecc(nu, e):
    """
    Convert true anomaly to eccentric anomaly.

    Parameters
    ----------
    nu : float
        True anomaly (radians).

    e : float
        Eccentricity.

    Returns
    -------
    E : float
        Eccentric anomaly (radians).
    """
    E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu / 2))
    return E


def true2mean(nu, e):
    """
    Convert true anomaly to mean anomaly.

    Parameters
    ----------
    nu : float
        True anomaly (radians).

    e : float
        Eccentricity.

    Returns
    -------
    M : float
        Mean anomaly (radians).
    """
    E = true2ecc(nu, e)
    M = ecc2mean(E, e)
    return M


def cart_to_coe(rv, GM):
    """
    Convert Cartesian Coordinates to Classical Orbital Elements

    Args:
        rv (np.array, Nx6): array of cartesian coordinates [x, y, z, vx, vy, vz]
    """
    r_vec = rv[0:3]
    v_vec = rv[3:6]

    # 1
    h_bar = np.cross(r_vec, v_vec)  # N x 3
    h = np.linalg.norm(h_bar)  #
    # 2
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    # 3
    E = 0.5 * (v**2) - GM / r
    # 4
    a = -GM / (2 * E)
    # 5
    e = np.sqrt(1 - (h**2) / (a * GM))
    # 6
    i = np.arccos(h_bar[2] / h)
    # 7
    omega_LAN = np.arctan2(h_bar[0], -h_bar[1])
    # 8
    # beware of division by zero here
    lat = np.arctan2(
        np.divide(r_vec[2], (np.sin(i))),
        (r_vec[0] * np.cos(omega_LAN) + r_vec[1] * np.sin(omega_LAN)),
    )
    # 9
    p = a * (1 - e**2)
    dot_rv = np.sum(r_vec * v_vec)
    nu = np.arctan2(np.sqrt(p / GM) * dot_rv, p - r)
    # 10
    omega_AP = lat - nu

    oe = np.zeros(6)
    oe[0] = a
    oe[1] = e
    oe[2] = i
    oe[3] = omega_LAN
    oe[4] = omega_AP
    EA = true2ecc(nu, e)
    MA = ecc2mean(EA, e)
    oe[5] = MA

    return oe


def coe_to_cart(x_oe, GM):
    """
    Convert classical orbital elements to Cartesian coordinates.

    Parameters
    ----------
    x_oe : array_like (6,) or (N, 6)
        Classical orbital elements [a, e, i, Omega, w, M]
    GM : float
        Gravitational parameter.

    Returns
    -------
    rv : array_like (6,) or (N, 6)
        Cartesian coordinates [x, y, z, vx, vy, vz]

    """
    a, e, i, Omega, w, M = x_oe.T

    # solve the kepler equation for E
    E = M + e * np.sin(M)  # Initial guess
    for _ in range(10):  # Iterate to refine the guess
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))

    # Create perifocal coordinate vectors
    P0 = np.cos(w) * np.cos(Omega) - np.sin(w) * np.cos(i) * np.sin(Omega)
    P1 = np.cos(w) * np.sin(Omega) + np.sin(w) * np.cos(i) * np.cos(Omega)
    P2 = np.sin(w) * np.sin(i)
    P = np.hstack((P0, P1, P2))

    Q0 = -np.sin(w) * np.cos(Omega) - np.cos(w) * np.cos(i) * np.sin(Omega)
    Q1 = -np.sin(w) * np.sin(Omega) + np.cos(w) * np.cos(i) * np.cos(Omega)
    Q2 = np.cos(w) * np.sin(i)
    Q = np.hstack((Q0, Q1, Q2))

    p = a * (1 - e * e)
    rnorm = a * (1 - e * np.cos(E))

    c1 = a * (np.cos(E) - e)
    c2 = np.sqrt(a * p) * np.sin(E)
    r = c1 * P + c2 * Q

    c3 = -(np.sqrt(GM * a) / rnorm * (np.sin(E)))
    c4 = np.sqrt(GM * p) / rnorm * np.cos(E)
    v = c3 * P + c4 * Q

    rv = np.hstack((r, v))

    return rv


def coe_to_mqoe(coe):
    """
    Convert classical orbital elements to modified equinoctial elements.

    Parameters
    ----------
    coe : array_like (6,) or (N, 6)
        Classical orbital elements [a, e, i, Omega, w, M].

    Returns
    -------
    mqoe : array_like (6,) or (N, 6)
        Modified equinoctial elements [a, h, k, p, q, M].
    """
    a, e, i, Omega, w, M = coe.T

    nu = mean_to_true(M, e)  # True anomaly

    p = a * (1 - e**2)
    f = e * np.cos(w + Omega)  # 0 -- 1
    g = e * np.sin(w + Omega)  # 0 -- 1
    h = np.tan(i / 2) * np.cos(Omega)
    k = np.tan(i / 2) * np.sin(Omega)
    L = wrapTo2Pi(Omega + w + nu)

    if coe.ndim == 1:
        return np.array([p, f, g, h, k, L])
    else:
        return np.vstack([p, f, g, h, k, L]).T


def mqoe_to_coe(mqoe, use_true_anom=False):
    """
    Convert modified equinoctial elements to classical orbital elements.

    Parameters
    ----------
    mqoe : array_like (6,) or (N, 6)
        Modified equinoctial elements [p, f, g, h, k, L].

    Returns
    -------
    coe : array_like (6,) or (N, 6)
        Classical orbital elements [a, e, i, Omega, w, M].
    """
    p, f, g, h, k, L = mqoe.T

    a = p / (1 - f**2 - g**2)
    e = np.sqrt(f**2 + g**2)
    i = np.arctan2(2 * np.sqrt(h**2 + k**2), 1 - h**2 - k**2)
    Omega = np.arctan2(k, h)
    w = np.arctan2(g, f) - Omega
    nu = L - w - Omega

    if use_true_anom:
        M = nu
    else:
        M = wrapTo2Pi(true2mean(nu, e))

    if mqoe.ndim == 1:
        coe = np.array([a, e, i, Omega, w, M])
    else:
        coe = np.vstack([a, e, i, Omega, w, M]).T

    return coe


def mqoe_to_cart(mqoe, GM):
    """
    Convert modified equinoctial elements to Cartesian coordinates.

    Parameters
    ----------
    mqoe : array_like (6,) or (N, 6)
        Modified equinoctial elements [p, f, g, h, k, L].

    GM : float
        Gravitational parameter.

    Returns
    -------
    rv : array_like (6,) or (N, 6)
        Cartesian coordinates [x, y, z, vx, vy, vz].
    """
    if mqoe.ndim == 1:
        mqoe = mqoe.reshape(1, -1)

    p, f, g, h, k, L = mqoe.T

    hh = h * h
    kk = k * k
    tkh = 2 * k * h
    s2 = 1 + hh + kk

    cL = np.cos(L)
    sL = np.sin(L)
    w = 1 + f * cL + g * sL
    r = p / w
    smp = np.sqrt(GM / p)

    # Build the two orthonormal frame vectors fhat, ghat
    fhat = np.vstack([1 - kk + hh, tkh, -2 * k]).T / s2[:, np.newaxis]  # N x 3
    ghat = np.vstack([tkh, 1 + kk - hh, 2 * h]).T / s2[:, np.newaxis]  # N x 3

    # In‐plane coordinates and rates
    x = r * cL
    y = r * sL
    xdot = -smp * (g + sL)
    ydot = smp * (f + cL)

    # Position and velocity in the inertial frame
    pos = x[:, np.newaxis] * fhat + y[:, np.newaxis] * ghat
    vel = xdot[:, np.newaxis] * fhat + ydot[:, np.newaxis] * ghat

    rv = np.hstack((pos, vel))

    if rv.shape[0] == 1:
        rv = rv[0]

    return rv


def cart_to_mqoe(rv, GM):
    """
    Convert Cartesian coordinates to modified equinoctial elements.

    Parameters
    ----------
    rv : array_like (6,) or (N, 6)
        Cartesian coordinates [x, y, z, vx, vy, vz].

    GM : float
        Gravitational parameter.

    Returns
    -------
    mqoe : array_like (6,) or (N, 6)
        Modified equinoctial elements [p, f, g, h, k, L].
    """
    # if rv is a 1D array, convert it to a 2D array
    if rv.ndim == 1:
        rv = rv.reshape(1, -1)

    # Split position & velocity
    r = rv[:, :3]
    v = rv[:, 3:6]

    # Norms and basic vectors
    rmag = np.linalg.norm(r, axis=1)
    rdv = np.sum(r * v, axis=1)  # dot product
    hvec = np.cross(r, v)
    hmag = np.linalg.norm(hvec, axis=1)

    # Unit vectors
    rhat = r / rmag[:, np.newaxis]
    hhat = hvec / hmag[:, np.newaxis]
    vhat = (rmag[:, np.newaxis] * v - rdv[:, np.newaxis] * rhat) / hmag[:, np.newaxis]

    # Semi-latus rectum
    p = hmag**2 / GM

    # Equinoctial inclination parameters
    # Note: in Fortran: k = hhat(1)/(1 + hhat(3)), h = -hhat(2)/(1 + hhat(3))
    k = hhat[:, 0] / (1.0 + hhat[:, 2])
    h = -hhat[:, 1] / (1.0 + hhat[:, 2])

    kk = k**2
    hh = h**2
    tkh = 2.0 * k * h
    s2 = 1.0 + hh + kk

    # Eccentricity vector
    ecc = np.cross(v, hvec) / GM - rhat

    # Frame vectors fhat, ghat
    fhat = np.vstack([1.0 - kk + hh, tkh, -2.0 * k]).T / s2[:, np.newaxis]

    ghat = np.vstack([tkh, 1.0 + kk - hh, 2.0 * h]).T / s2[:, np.newaxis]

    # Equinoctial shape parameters
    f = np.sum(ecc * fhat, axis=1)  #
    g = np.sum(ecc * ghat, axis=1)

    # True longitude L
    # Fortran: L = atan2(rhat(2)-vhat(1), rhat(1)+vhat(2))
    L = np.arctan2(rhat[:, 1] - vhat[:, 0], rhat[:, 0] + vhat[:, 1])
    L = wrapTo2Pi(L)

    meoe = np.vstack([p, f, g, h, k, L]).T

    if meoe.shape[0] == 1:
        meoe = meoe[0]

    return meoe


def rot_moon_ci2pa(t_tai, rot_only=False):

    angles = pnt.get_lunar_orientation_angles(t_tai)  # MCI to PA angles
    phi, theta, psi, phi_dot, theta_dot, psi_dot = angles

    R_z_phi = RotZ(phi)
    R_x_theta = RotX(theta)
    R_z_psi = RotZ(psi)

    R_mi2pa = R_z_psi @ R_x_theta @ R_z_phi

    if rot_only:
        rot_mat = np.zeros((6, 6))
        rot_mat[0:3, 0:3] = R_mi2pa
        rot_mat[3:6, 3:6] = R_mi2pa
    else:
        R_mi2pa_dot = (
            RotZ_dot(psi, psi_dot) @ R_x_theta @ R_z_phi
            + R_z_psi @ RotX_dot(theta, theta_dot) @ R_z_phi
            + R_z_psi @ R_x_theta @ RotZ_dot(phi, phi_dot)
        )
        rot_mat = np.zeros((6, 6))
        rot_mat[0:3, 0:3] = R_mi2pa
        rot_mat[3:6, 3:6] = R_mi2pa
        rot_mat[3:6, 0:3] = R_mi2pa_dot

    return rot_mat


def rot_moon_pa2ci(t_tai, rot_only=False):
    """
    Inverse of rot_moon_ci2pa: maps PA state -> CI state.
    Returns a 6x6 matrix with position and velocity blocks.
    """
    # Same angles (they define CI->PA), we just transpose to invert
    phi, theta, psi, phi_dot, theta_dot, psi_dot = pnt.get_lunar_orientation_angles(t_tai)

    R_z_phi = RotZ(phi)
    R_x_theta = RotX(theta)
    R_z_psi = RotZ(psi)

    # CI -> PA rotation
    R_mi2pa = R_z_psi @ R_x_theta @ R_z_phi
    # PA -> CI rotation is the transpose
    R_pa2mi = R_mi2pa.T

    if rot_only:
        rot_mat = np.zeros((6, 6))
        rot_mat[0:3, 0:3] = R_pa2mi
        rot_mat[3:6, 3:6] = R_pa2mi
        return rot_mat

    # Time derivative of CI->PA rotation (already in your forward function)
    R_mi2pa_dot = (
        RotZ_dot(psi, psi_dot) @ R_x_theta @ R_z_phi
        + R_z_psi @ RotX_dot(theta, theta_dot) @ R_z_phi
        + R_z_psi @ R_x_theta @ RotZ_dot(phi, phi_dot)
    )

    # Cross term for the inverse:
    # d/dt(R^{-1}) = -R^{-1} R_dot R^{-1}  ⇒  (R_pa2mi)_dot = -R^T R_dot R^T
    cross = -R_pa2mi @ R_mi2pa_dot @ R_pa2mi

    rot_mat = np.zeros((6, 6))
    rot_mat[0:3, 0:3] = R_pa2mi  # position block
    rot_mat[3:6, 3:6] = R_pa2mi  # velocity rotation
    rot_mat[3:6, 0:3] = cross  # velocity-position coupling

    return rot_mat


def convert_ci2pa(t_tai_vec, rv_vec, rotate_only=False):
    """
    Convert from CI to PA frame.

    Parameters
    ----------
    t_tai_vec : array_like
        Time vector in TAI.
    rv_vec : array_like
        Position and velocity vector in CI frame.
    rot_only : bool, optional
        If True, return only the rotation matrix. Default is False.

    Returns
    -------
    rv_pa_vec : array_like
        Position and velocity vector in PA frame.
    """
    if rv_vec.ndim == 1:
        rv_vec = rv_vec.reshape(1, -1)
    if np.isscalar(t_tai_vec):
        t_tai_vec = np.array([t_tai_vec])

    rv_pa_vec = np.zeros_like(rv_vec)

    for i, t_tai in enumerate(t_tai_vec):
        R_mi2pa = rot_moon_ci2pa(t_tai, rot_only=rotate_only)
        rv_pa_vec[i] = R_mi2pa @ rv_vec[i]

    if rv_pa_vec.shape[0] == 1:
        rv_pa_vec = rv_pa_vec[0]  # return 1D array if input was 1D

    return rv_pa_vec


def convert_pa2ci(t_tai_vec, rv_vec, rotate_only=False):
    """
    Convert from PA to CI frame.

    Parameters
    ----------
    t_tai_vec : array_like
        Time vector in TAI.
    rv_vec : array_like
        Position and velocity vector in PA frame.
    rot_only : bool, optional
        If True, return only the rotation matrix. Default is False.

    Returns
    -------
    rv_ci_vec : array_like
        Position and velocity vector in CI frame.
    """
    if rv_vec.ndim == 1:
        rv_vec = rv_vec.reshape(1, -1)
    if np.isscalar(t_tai_vec):
        t_tai_vec = np.array([t_tai_vec])

    rv_ci_vec = np.zeros_like(rv_vec)

    for i, t_tai in enumerate(t_tai_vec):
        R_pa2mi = rot_moon_pa2ci(t_tai, rot_only=rotate_only)
        rv_ci_vec[i] = R_pa2mi @ rv_vec[i]

    if rv_ci_vec.shape[0] == 1:
        rv_ci_vec = rv_ci_vec[0]  # return 1D array if input was 1D

    return rv_ci_vec
