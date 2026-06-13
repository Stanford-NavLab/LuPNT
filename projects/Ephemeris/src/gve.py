import pylupnt as pnt
import numpy as np
from scipy.integrate import solve_ivp
from copy import deepcopy

try:
    import src.gravity_field as gf
except ImportError:
    import gravity_field as gf


def get_rtn_matrix(rv):
    r = rv[:3]
    v = rv[3:]
    r_hat = r / np.linalg.norm(r)
    h = np.cross(r, v)
    h_hat = h / np.linalg.norm(h)
    t_hat = np.cross(h_hat, r_hat)
    t_hat = t_hat / np.linalg.norm(t_hat)
    M = np.vstack([r_hat, t_hat, h_hat])
    return M


def gauss_variational_equation(x_oe, a_rtn, mu):
    """
    Compute the time derivatives of the osculating orbital elements using
    Gauss' variational equations.

    Reference:
    Battin "An Introduction to the Mathematics and Methods of Astrodynamics",
    equation (10.41), page 488

    Parameters:
    -----------
    x_oe : array_like, shape (6,)
        Osculating orbital elements [a, e, i, Omega, omega, f] where:
          a     = semi-major axis
          e     = eccentricity
          i     = inclination (radians)
          Omega = RAAN (radians)
          omega = argument of perigee (radians)
          f   = true anomaly (radians)
    a_rtn : array_like, shape (3,)
        Acceleration in the RTN frame [u_r, u_t, u_h].
    mu : float
        Gravitational parameter of the central body.

    Returns:
    --------
    dx_dt : ndarray, shape (6,)
        Time derivatives of the osculating elements.
    """
    # Unpack orbital elements
    a, e, inc, Omega, omega, f = x_oe

    # Compute auxiliary quantities
    theta = omega + f  # Argument relative to node line in orbital plane
    h = np.sqrt(mu * a * (1 - e**2))  # Specific angular momentum
    p = a * (1 - e**2)  # Semi-latus rectum
    r = p / (1 + e * np.cos(f))  # Radial distance
    # v = np.sqrt(mu * (2 / r - 1 / a))  # Orbital velocity

    # convert from rtn to polar frame
    u_r, u_th, u_h = a_rtn

    # Gauss variational equations

    # Semi-major axis derivative
    da_dt = (2 * a**2 / h) * (e * np.sin(f) * u_r + (p / r) * u_th)

    # Eccentricity derivative
    de_dt = 1 / h * (p * np.sin(f) * u_r + ((p + r) * np.cos(f) + r * e) * u_th)

    # Inclination derivative
    di_dt = (r * np.cos(theta) / h) * u_h

    # Right Ascension of the Ascending Node derivative
    dOmega_dt = (r * np.sin(theta) / (h * np.sin(inc))) * u_h

    # Argument of perigee derivative
    domega_dt = (1 / (h * e)) * (-p * np.cos(f) * u_r + (p + r) * np.sin(f) * u_th) - np.cos(
        inc
    ) * dOmega_dt

    # True anomaly derivative
    df_dt = (h / r**2) + (1 / (h * e)) * (p * np.cos(f) * u_r - (p + r) * np.sin(f) * u_th)

    # Pack derivatives into a single array
    dx_dt = np.array([da_dt, de_dt, di_dt, dOmega_dt, domega_dt, df_dt])
    return dx_dt


###############################################################################
# Gauss Variational Equations
###############################################################################
class GVE:
    def __init__(self, third_bodies=[pnt.EARTH], sph=[2, 2], center=pnt.MOON):
        self.third_bodies = third_bodies

        bodies_GM = []
        for body in self.third_bodies:
            if body == pnt.EARTH:
                bodies_GM.append(pnt.GM_EARTH)
            elif body == pnt.SUN:
                bodies_GM.append(pnt.GM_SUN)
            elif body == pnt.MOON:
                bodies_GM.append(pnt.GM_MOON)
            else:
                # error
                print("Unsupported third body")
                return

        self.bodies_GM = bodies_GM

        self.sph = sph  # Spherical harmonics model (order, degree)
        if len(sph) != 2:
            # error
            print("Spherical harmonics model must be a 2-element list")
            return

        self.center = center
        if center == pnt.MOON:
            self.mu = pnt.GM_MOON
            self.frame_i = pnt.MOON_CI
            self.frame_b = pnt.MOON_PA
            self.gravfile = "grgm900c.cof"
        else:
            # error
            print("Unsupported central body: only Moon is supported")
            return

        self.grav = gf.read_harmonic_gravity_field(
            pnt.find_file(self.gravfile), n=sph[0], m=sph[1], normalized=True
        )

    def compute_frtn(self, t_tai, x_oe):
        """
        Compute the acceleration in the RTN frame using Gauss' variational equations.

        Parameters:
        -----------
        t_tai : float
            Time in TAI (seconds).
        x_oe : array_like
            Osculating elements [a, e, i, Omega, omega, f] where:
              a     = semi-major axis
              e     = eccentricity
              i     = inclination (radians)
              Omega = RAAN (radians)
              omega = argument of perigee (radians)
              f   = true anomaly (radians)
        Returns:
        --------
        acc_rtn : ndarray, shape (3,)
            Acceleration in the RTN frame [u_r, u_t, u_h].

        """

        x_oe_M = deepcopy(x_oe)
        # x_oe_M[5] = true2mean(x_oe[5], x_oe[1]) # convert to true anomaly
        # x_cart = coe_to_cart(x_oe_M, self.mu)
        x_oe_M[5] = pnt.true2mean_anomaly(x_oe[5], x_oe[1])  # convert to true anomaly
        x_cart = pnt.classical_to_cart(x_oe_M, self.mu)

        # third body
        acc_body = np.zeros(3)
        if len(self.third_bodies) == 0:
            acc_body = np.zeros(3)
        else:
            for i, body in enumerate(self.third_bodies):
                rv_body = pnt.get_body_pos_vel(t_tai, body, self.frame_i)
                acc_body += pnt.acceleration_point_mass(
                    x_cart[:3], rv_body[:3], GM=self.bodies_GM[i]
                )

        # high-order gravity field
        if (self.sph[0] > 0) or (self.sph[1] > 0):
            r_bf = pnt.convert_frame(t_tai, x_cart[:3], self.frame_i, self.frame_b)
            grav = self.grav
            CS = grav.CS.copy()  # Copy the coefficients
            CS[0, 0] = 0.0  # remove the monopole term
            acc_grav_bf = pnt.acceleration_gravity_field(r_bf, grav.GM, grav.R, CS, grav.n, grav.m)
            # subtract the central body acceleration
            # do not include the acceleration due to the central body
            acc_grav_i = pnt.convert_frame(t_tai, acc_grav_bf, self.frame_b, self.frame_i)
        else:
            acc_grav_i = np.zeros(3)

        # total acceleration
        acc_i = acc_grav_i + acc_body  # total acceleration in the inertial frame
        acc_rtn = get_rtn_matrix(x_cart) @ acc_i  # total acceleration in the RTN frame

        return acc_rtn

    def propagate(self, t_span, x_oe0, method="RK45"):
        """
        Propagate the osculating elements using Gauss variational equations.

        Parameters:
        -----------
        t_span : array_like
            Time span for the propagation.
        x_oe0 : array_like
            Initial osculating elements. (a, e, i, Omega, omega, M: mean anomaly)
        method : str
            Integration method. Default is 'RK45'.

        Returns:
        --------
        t : array_like
            Time vector.
        x_oe : array_like
            Osculating elements at each time step.
        """

        t0 = t_span[0]
        tf = t_span[-1]
        tprop = [t0, tf]

        # first convert x_oe0 to true anomaly (GVE is propagated in true anomaly)
        x_oe0 = deepcopy(x_oe0)
        x_oe0[5] = pnt.mean_to_true_anomaly(x_oe0[5], x_oe0[1])

        # ODE function
        def f_ode(t, x):
            # x: osculating elements
            a_rtn = self.compute_frtn(t, x)
            dx = gauss_variational_equation(x, a_rtn, mu=self.mu)
            return dx

        # integrate the ODE
        sol = solve_ivp(f_ode, tprop, x_oe0, method=method, t_eval=t_span, rtol=1e-13, atol=1e-13)

        t = sol.t  # time points
        x_oe = sol.y.T  # Transpose to have the shape (n, 6)

        # convert to mean anomaly
        for i, coe in enumerate(x_oe):
            x_oe[i, 5] = pnt.true2mean_anomaly(x_oe[i, 5], x_oe[i, 1])

        return t, x_oe
