import numpy as np
import pylupnt as pnt
from scipy.integrate import solve_ivp

try:
    import src.gravity_field as gf
except ImportError:
    import gravity_field as gf


class NBodyProp:
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

    def compute_f(self, t_tai, x_cart):

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
            acc_grav_bf = pnt.acceleration_gravity_field(
                r_bf, grav.GM, grav.R, grav.CS, grav.n, grav.m
            )
            acc_grav_i = pnt.convert_frame(t_tai, acc_grav_bf, self.frame_b, self.frame_i)
        else:
            acc_grav_i = (
                -self.mu * x_cart[:3] / np.linalg.norm(x_cart[:3]) ** 3
            )  # central body acceleration

        # total acceleration
        acc_i = acc_grav_i + acc_body  # total acceleration in the inertial frame

        return acc_i

    def propagate(self, t_span, x_cart0, method="RK45"):
        """
        Propagate the osculating elements using Gauss variational equations.

        Parameters:
        -----------
        t_span : array_like
            Time span for the propagation.
        x_cart0 : array_like
            Initial cartesian state vector [r, v].
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

        def f_ode(t, rv):
            # x: osculating elements
            acc = self.compute_f(t, rv)
            v = rv[3:]
            dx = np.concatenate([v, acc])
            return dx

        sol = solve_ivp(f_ode, tprop, x_cart0, method=method, t_eval=t_span, rtol=1e-13, atol=1e-13)

        t = sol.t  # time points
        x_oe = sol.y.T  # Transpose to have the shape (n, 6)

        return t, x_oe
