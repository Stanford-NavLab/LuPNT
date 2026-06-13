import numpy as np
import pylupnt as pnt
import matplotlib.pyplot as plt
import copy
import cvxpy as cp
from scipy.linalg import lstsq

try:
    from .ephemeris import Ephemeris
    from .basis import get_poly_basis, get_poly_basis_dt, get_fourier_basis, get_fourier_basis_dt
    from .orbit_utils import wrapToPi, wrapTo2Pi, rms
except ImportError:
    from ephemeris import Ephemeris
    from basis import get_poly_basis, get_poly_basis_dt, get_fourier_basis, get_fourier_basis_dt
    from orbit_utils import wrapToPi, wrapTo2Pi, rms


class KeplarianEphemeris(Ephemeris):
    def __init__(self, param_dict, body=pnt.MOON, print_info=False):
        super().__init__(body)
        self.param_dict = param_dict

        self.create_ephem_dict(print_info=print_info)

    def create_ephem_dict(self, print_info=False):

        self.params = ["a", "e", "w", "M", "f", "r", "u", "i", "l"]
        self.idx_keys = {}
        self.idx_keys["t_ref"] = 0
        self.idx_keys["t_fit"] = 1

        self.keys_list = []
        self.keys_list.append("t_ref")
        self.keys_list.append("t_fit")

        self.is_const = {}
        self.poly_orders = {}
        self.fourier_orders = {}
        self.fourier_sidereal_orders = {}
        self.poly_types = {}
        self.use_poly = False  # any parameter use polynomial
        self.use_cheby = False  # any parameter use chebyshev
        self.use_legendre = False  # any parameter use legendre
        self.use_fourier = False  # any parameter use fourier
        self.use_fourier_sidereal = False  # any parameter use fourier sidereal
        self.max_poly_order = 0
        self.max_fourier_order = 0
        self.max_fourier_sidereal_order = 0

        for param_str in self.params:
            self.is_const[param_str] = False

            if param_str in self.param_dict.keys():
                param = self.param_dict[param_str]
                self.poly_orders[param_str] = param["linear"]
                self.fourier_orders[param_str] = param["fourier"]
                self.fourier_sidereal_orders[param_str] = param.get("fourier_sidereal", 0)
                self.poly_types[param_str] = param["polytype"]

                self.max_poly_order = max(self.max_poly_order, self.poly_orders[param_str])
                self.max_fourier_order = max(self.max_fourier_order, self.fourier_orders[param_str])
                self.max_fourier_sidereal_order = max(
                    self.max_fourier_sidereal_order, self.fourier_sidereal_orders[param_str]
                )

                if self.poly_types[param_str] == "monomial":
                    self.use_poly = True
                if self.poly_types[param_str] == "chebyshev":
                    self.use_cheby = True
                if self.poly_types[param_str] == "legendre":
                    self.use_legendre = True

                # linear orders
                for i in range(self.poly_orders[param_str] + 1):
                    key = f"{param_str}_{i}"
                    self.idx_keys[key] = len(self.keys_list)
                    self.keys_list.append(key)

                # fourier orders
                if self.fourier_orders[param_str] > 0:
                    self.use_fourier = True
                    for i in range(self.fourier_orders[param_str]):
                        key = f"{param_str}_C{i+1}"
                        self.idx_keys[key] = len(self.keys_list)
                        self.keys_list.append(key)
                        key = f"{param_str}_S{i+1}"
                        self.idx_keys[key] = len(self.keys_list)
                        self.keys_list.append(key)

                # fourier sidereal orders
                if self.fourier_sidereal_orders[param_str] > 0:
                    self.use_fourier_sidereal = True
                    for i in range(self.fourier_sidereal_orders[param_str]):
                        key = f"{param_str}_Cs{i+1}"
                        self.idx_keys[key] = len(self.keys_list)
                        self.keys_list.append(key)
                        key = f"{param_str}_Ss{i+1}"
                        self.idx_keys[key] = len(self.keys_list)
                        self.keys_list.append(key)

                # check if the parameter is constant
                if (
                    self.poly_orders[param_str] == 0
                    and self.fourier_orders[param_str] == 0
                    and self.fourier_sidereal_orders[param_str] == 0
                ):
                    self.is_const[param_str] = True

        self.n_params = len(self.keys_list)

        if print_info:
            print("key list: ", self.keys_list)
            print("keys idx: ", self.idx_keys)
            print("Number of parameters: ", self.n_params)
            for param_str in self.params:
                print(
                    f"{param_str}: linear order: {self.poly_orders[param_str]}, fourier order: {self.fourier_orders[param_str]}, fourier sidereal order: {self.fourier_sidereal_orders[param_str]}, poly type: {self.poly_types[param_str]}"
                )

        return

    def coe2dict(self, coe, t_ref, t_fit):
        dict = {}

        for i, key in enumerate(self.keys_list):
            if key == "t_ref":
                dict[key] = t_ref
            elif key == "t_fit":
                dict[key] = t_fit
            elif key == "a_0":
                dict[key] = coe[0]
            elif key == "e_0":
                dict[key] = coe[1]
            elif key == "i_0":
                dict[key] = coe[2]
            elif key == "l_0":
                dict[key] = coe[3]
            elif key == "w_0":
                dict[key] = coe[4]
            elif key == "M_0":
                dict[key] = coe[5]
            else:
                dict[key] = 0

        return dict

    def print_ephem(self, ephem):
        ephem_dict = self.ephem2dict(ephem)
        for key, value in ephem_dict.items():
            if abs(value) >= 0.01:
                print("{0:}: {1:.4f}".format(key, value))
            else:
                print("{0:}: {1:.2e}".format(key, value))

    def get_value(self, ephem_dict, key, default=0):
        """
        Get the value of the key in the ephemeris
        """
        if key in ephem_dict.keys():
            return ephem_dict[key]
        else:
            return default

    def get_poly_coeffs(self, dict, param_str):
        """
        Get the coefficients of the polynomial
        """
        coeffs = []
        for i in range(self.poly_orders[param_str] + 1):
            param_str_i = f"{param_str}_{i}"
            coeffs.append(self.get_value(dict, param_str_i, 0))
        return np.array(coeffs)

    def get_poly_coeff_keys(self, param_str):
        """
        Get the keys of the polynomial
        """
        keys = []
        for i in range(self.poly_orders[param_str] + 1):
            param_str_i = f"{param_str}_{i}"
            keys.append(param_str_i)
        return keys

    def get_fourier_coeffs(self, dict, param_str):
        """
        Get the coefficients of the fourier
        """
        coeffs = []
        for i in range(self.fourier_orders[param_str]):
            param_str_C = f"{param_str}_C{i+1}"
            param_str_S = f"{param_str}_S{i+1}"
            coeffs.append(self.get_value(dict, param_str_C, 0))
            coeffs.append(self.get_value(dict, param_str_S, 0))
        return np.array(coeffs)

    def get_fourier_coeff_keys(self, param_str):
        """
        Get the keys of the fourier
        """
        keys = []
        for i in range(self.fourier_orders[param_str]):
            param_str_C = f"{param_str}_C{i+1}"
            param_str_S = f"{param_str}_S{i+1}"
            keys.append(param_str_C)
            keys.append(param_str_S)
        return keys

    def get_fourier_sidereal_coeffs(self, dict, param_str):
        """
        Get the coefficients of the fourier sidereal
        """
        coeffs = []
        for i in range(self.fourier_sidereal_orders[param_str]):
            param_str_Cs = f"{param_str}_Cs{i+1}"
            param_str_Ss = f"{param_str}_Ss{i+1}"
            coeffs.append(self.get_value(dict, param_str_Cs, 0))
            coeffs.append(self.get_value(dict, param_str_Ss, 0))
        return np.array(coeffs)

    def get_fourier_sidereal_coeff_keys(self, param_str):
        """
        Get the keys of the fourier sidereal
        """
        keys = []
        for i in range(self.fourier_sidereal_orders[param_str]):
            param_str_Cs = f"{param_str}_Cs{i+1}"
            param_str_Ss = f"{param_str}_Ss{i+1}"
            keys.append(param_str_Cs)
            keys.append(param_str_Ss)
        return keys

    def compute_param(self, dict, param_str, t_k, Phi, offset):
        """
        Compute the value of the parameter at time t_k
        """

        pvalue = offset

        # sidereal period
        omega_sidereal = 2 * np.pi / self.T_sidereal

        # pre-computation for chebyshev
        if self.use_poly:
            Tp = get_poly_basis(t_k, self.max_poly_order, dict["t_fit"], "monomial")
        if self.use_cheby:
            Tc = get_poly_basis(t_k, self.max_poly_order, dict["t_fit"], "chebyshev")
        if self.use_fourier and Phi is not None:
            Tf = get_fourier_basis(Phi, self.max_fourier_order)
        if self.use_fourier_sidereal:
            Tfs = get_fourier_basis(omega_sidereal * t_k, self.max_fourier_sidereal_order)

        if param_str in self.params:
            ln = self.poly_orders[param_str]
            fn = self.fourier_orders[param_str]
            fsn = self.fourier_sidereal_orders[param_str]
            ptype = self.poly_types[param_str]

            # linear terms
            l_coeffs = self.get_poly_coeffs(dict, param_str)
            if ptype == "monomial":
                pvalue += Tp[:, : (ln + 1)] @ l_coeffs
            elif ptype == "chebyshev":
                pvalue += Tc[:, : (ln + 1)] @ l_coeffs
            else:
                print("Invalid poly type: ", ptype)
                return None

            # fourier terms
            if fn > 0:
                f_coeffs = self.get_fourier_coeffs(dict, param_str)
                pvalue += Tf[:, : (fn * 2)] @ f_coeffs

            # fourier sidereal terms
            if fsn > 0:
                fs_coeffs = self.get_fourier_sidereal_coeffs(dict, param_str)
                pvalue += Tfs[:, : (fsn * 2)] @ fs_coeffs

        else:
            print("Invalid key: ", param_str)
            return None

        return pvalue

    def compute_param_dt(self, dict, param_str, t_k, Phi, offset):
        """
        Compute the value of the parameter at time t_k
        """

        pvalue = offset

        # sidereal period
        omega_sidereal = 2 * np.pi / self.T_sidereal

        # pre-computation for chebyshev
        if self.use_poly:
            Tpdot = get_poly_basis_dt(t_k, self.max_poly_order, dict["t_fit"], "monomial")
        if self.use_cheby:
            Tcdot = get_poly_basis_dt(t_k, self.max_poly_order, dict["t_fit"], "chebyshev")
        if self.use_fourier and Phi is not None:
            Phidot = self.compute_Phidot(dict, t_k)
            Tfdot = get_fourier_basis_dt(Phi, Phidot, self.max_fourier_order)
        if self.use_fourier_sidereal:
            Tfsdot = get_fourier_basis_dt(
                omega_sidereal * t_k, omega_sidereal, self.max_fourier_sidereal_order
            )

        if param_str in self.params:
            ln = self.poly_orders[param_str]
            fn = self.fourier_orders[param_str]
            fsn = self.fourier_sidereal_orders[param_str]
            ptype = self.poly_types[param_str]

            # linear terms
            l_coeffs = self.get_poly_coeffs(dict, param_str)
            if ptype == "monomial":
                pvalue += Tpdot[:, : (ln + 1)] @ l_coeffs
            elif ptype == "chebyshev":
                pvalue += Tcdot[:, : (ln + 1)] @ l_coeffs
            else:
                print("Invalid poly type: ", ptype)
                return None

            # fourier terms
            if fn > 0 and Phi is not None:
                f_coeffs = self.get_fourier_coeffs(dict, param_str)
                pvalue += Tfdot[:, : (fn * 2)] @ f_coeffs

            # fourier sidereal terms
            if fsn > 0:
                fs_coeffs = self.get_fourier_sidereal_coeffs(dict, param_str)
                pvalue += Tfsdot[:, : (fsn * 2)] @ fs_coeffs

        else:
            print("Invalid key: ", param_str)
            return None

        return pvalue

    def compute_Phidot(self, dict, t_k):

        a_0 = dict["a_0"]
        T_orbit = 2 * np.pi * np.sqrt(a_0**3 / self.GM)
        Phi_tmp = 2 * np.pi / T_orbit * t_k  # tempoaral variable for fourier

        a_k = self.compute_param(dict, "a", t_k, Phi_tmp, 0)
        n = np.sqrt(self.GM / a_k**3)  # Todo: replace this with integration (trapz?)
        Mdot = self.compute_param_dt(dict, "M", t_k, None, n)  # None to avoid circular dependency
        e = self.compute_param(dict, "e", t_k, Phi_tmp, 0)
        edot = self.compute_param_dt(dict, "e", t_k, None, 0)  # None to avoid circular dependency
        M = self.compute_param(dict, "M", t_k, Phi_tmp, n * t_k)

        E_k = M
        for _ in range(5):
            E_k = E_k - (E_k - e * np.sin(E_k) - M) / (1 - e * np.cos(E_k))

        Edot = (Mdot + edot * np.sin(E_k)) / (1 - e * np.cos(E_k))

        f = np.sqrt(1 - e * e) * np.sin(E_k)
        fdot = -e * edot / np.sqrt(1 - e * e) * np.sin(E_k) + np.sqrt(1 - e * e) * Edot * np.cos(
            E_k
        )
        g = 1 - e * np.cos(E_k)
        gdot = -edot * np.cos(E_k) + e * Edot * np.sin(E_k)
        cosnu = (np.cos(E_k) - e) / g
        nudot = (fdot * g - f * gdot) / (g * g) / cosnu

        wdot = self.compute_param_dt(dict, "w", t_k, None, 0)
        Phidot = nudot + wdot

        return Phidot

    def compute_rdot(self, dict, t_k):
        """
        Compute the radius at time t_k
        """
        a_0 = dict["a_0"]
        T_orbit = 2 * np.pi * np.sqrt(a_0**3 / self.GM)
        Phi_tmp = 2 * np.pi / T_orbit * t_k  # tempoaral variable for fourier

        a_k = self.compute_param(dict, "a", t_k, Phi_tmp, 0)
        n = np.sqrt(
            self.GM * np.ones_like(a_k) / a_k**3
        )  # Todo: replace this with integration (trapz?)
        e = self.compute_param(dict, "e", t_k, Phi_tmp, 0)
        M = self.compute_param(dict, "M", t_k, Phi_tmp, n * t_k)
        adot = self.compute_param_dt(dict, "a", t_k, Phi_tmp, 0)
        edot = self.compute_param_dt(dict, "e", t_k, Phi_tmp, 0)

        E_k = M
        for _ in range(5):
            E_k = E_k - (E_k - e * np.sin(E_k) - M) / (1 - e * np.cos(E_k))

        Mdot = self.compute_param_dt(dict, "M", t_k, Phi_tmp, n)

        # M = E - e sin(E)
        # Mdot = Edot - e * np.cos(E_k) * Edot - edot * np.sin(E_k)
        Edot = (Mdot + edot * np.sin(E_k)) / (1 - e * np.cos(E_k))

        # r = a * (1 - e * np.cos(E_k))
        rdot = adot - (
            edot * a_k * np.cos(E_k) + e * adot * np.cos(E_k) - e * a_k * Edot * np.sin(E_k)
        )

        return rdot

    def ephem2ruil(self, t, ephem, compute_velocity=False):
        dict = self.ephem2dict(ephem)

        t_ref = dict["t_ref"]
        t_k = t - t_ref

        a_0 = dict["a_0"]
        T_orbit = 2 * np.pi * np.sqrt(a_0**3 / self.GM)
        Phi_tmp = 2 * np.pi / T_orbit * t_k  # tempoaral variable for fourier

        a_k = self.compute_param(dict, "a", t_k, Phi_tmp, 0)

        n = np.sqrt(
            self.GM * np.ones_like(a_k) / a_k**3
        )  # Todo: replace this with integration (trapz?)
        M_k = self.compute_param(dict, "M", t_k, Phi_tmp, n * t_k)

        e_k = self.compute_param(dict, "e", t_k, Phi_tmp, 0)  # eccentricity
        # solve Kepler's equation E = M + e * sin(E)
        E_k = M_k
        for _ in range(5):
            E_k = E_k - (E_k - e_k * np.sin(E_k) - M_k) / (1 - e_k * np.cos(E_k))
        nu_k = 2 * np.arctan(np.sqrt((1 + e_k) / (1 - e_k)) * np.tan(E_k / 2))  # true anomaly

        # correct nu_k
        nu_k = self.compute_param(dict, "f", t_k, Phi_tmp, nu_k)

        # compute the argument of latitude
        w_k = self.compute_param(dict, "w", t_k, Phi_tmp, 0)  # argument of perigee
        Phi = w_k + nu_k

        # argument of latitude
        offset_u = w_k + nu_k
        u_k = self.compute_param(dict, "u", t_k, Phi, offset_u)  # argument of latitude

        # radius
        offset_r = a_k * (1 - e_k * np.cos(E_k))
        r_k = self.compute_param(dict, "r", t_k, Phi, offset_r)  # radius

        # inclination
        i_k = self.compute_param(dict, "i", t_k, Phi, 0)  # inclination

        # longitude of ascending node
        offset_l = -self.omega_b * t_k
        lambda_k = self.compute_param(dict, "l", t_k, Phi, offset_l)  # longitude of ascending node
        ruil = np.vstack([r_k, u_k, i_k, lambda_k]).T
        if compute_velocity:
            # for r and u, use diff
            rdot0 = self.compute_rdot(dict, t_k)
            Phidot = self.compute_Phidot(dict, t_k)
            drdt = self.compute_param_dt(dict, "r", t_k, Phi, rdot0)  # radius
            dudt = self.compute_param_dt(dict, "u", t_k, Phi, Phidot)  # argument of latitude
            didt = self.compute_param_dt(dict, "i", t_k, Phi, 0)  # inclination
            dldt = self.compute_param_dt(
                dict, "l", t_k, Phi, -self.omega_b
            )  # longitude of ascending node
            ruil_dot = np.vstack([drdt, dudt, didt, dldt]).T
        else:
            ruil_dot = np.zeros((1, 4))
        return ruil, ruil_dot

    def ephem2cart(self, t, ephem, compute_velocity=False, return_params=False):
        ruil, ruil_dot = self.ephem2ruil(t, ephem, compute_velocity)

        # if dimension is 1, we need to reshape it
        if ruil.ndim == 1:
            ruil = ruil.reshape(1, -1)
            ruil_dot = ruil_dot.reshape(1, -1)

        r_k = ruil[:, 0]
        u_k = ruil[:, 1]
        i_k = ruil[:, 2]
        lambda_k = ruil[:, 3]

        # positions
        x_k = r_k * np.cos(u_k)
        y_k = r_k * np.sin(u_k)
        x_pos = x_k * np.cos(lambda_k) - y_k * np.cos(i_k) * np.sin(lambda_k)
        y_pos = x_k * np.sin(lambda_k) + y_k * np.cos(i_k) * np.cos(lambda_k)
        z_pos = y_k * np.sin(i_k)

        # velocities
        if compute_velocity:
            rdot = ruil_dot[:, 0]
            udot = ruil_dot[:, 1]
            didt = ruil_dot[:, 2]
            lambdadot = ruil_dot[:, 3]

            xdot_k = rdot * np.cos(u_k) - r_k * udot * np.sin(u_k)
            ydot_k = rdot * np.sin(u_k) + r_k * udot * np.cos(u_k)

            x_vel = (
                -x_k * lambdadot * np.sin(lambda_k)
                - y_k
                * (
                    lambdadot * np.cos(i_k) * np.cos(lambda_k)
                    - didt * np.sin(i_k) * np.sin(lambda_k)
                )
                + xdot_k * np.cos(lambda_k)
                - ydot_k * np.cos(i_k) * np.sin(lambda_k)
            )
            y_vel = (
                x_k * lambdadot * np.cos(lambda_k)
                - y_k
                * (
                    lambdadot * np.cos(i_k) * np.sin(lambda_k)
                    + didt * np.sin(i_k) * np.cos(lambda_k)
                )
                + xdot_k * np.sin(lambda_k)
                + ydot_k * np.cos(i_k) * np.cos(lambda_k)
            )
            z_vel = y_k * didt * np.cos(i_k) + ydot_k * np.sin(i_k)

        # outputs
        res = np.vstack([x_pos, y_pos, z_pos]).T  # T x 3

        if compute_velocity:
            res = np.vstack([x_pos, y_pos, z_pos, x_vel, y_vel, z_vel]).T

        return res

    def fit_coeffs(
        self,
        residuals,
        param_str,
        ephem,
        t_k,
        Phi,
        residuals_dt=None,
        fix_endpoints=False,
        print_result=False,
    ):
        """
        Fit the coefficients of the parameter
        """
        n = t_k.size
        dict = self.ephem2dict(ephem)

        if self.is_const[param_str]:
            # if the parameter is constant, we don't need to fit it
            if print_result:
                print("  [{} (constant)] ".format(param_str))
                coeff = ephem[self.get_index(param_str + "_0")]
                res_rms = rms(residuals - coeff)
                if param_str == "M":
                    res_rms = np.mod(res_rms, np.pi * 2)
                print("    residual: ", res_rms)
                print("    coeffs  : ", coeff)
            return ephem

        # linear coefficients
        Amat = get_poly_basis(
            t_k, self.poly_orders[param_str], dict["t_fit"], self.poly_types[param_str]
        )
        keys = self.get_poly_coeff_keys(param_str)

        if self.fourier_orders[param_str] > 0:
            Amat_p = get_fourier_basis(Phi, self.fourier_orders[param_str])
            Amat = np.hstack((Amat, Amat_p))
            keys.extend(self.get_fourier_coeff_keys(param_str))

        if self.fourier_sidereal_orders[param_str] > 0:
            # sidereal period
            omega_sidereal = 2 * np.pi / self.T_sidereal
            Amat_ps = get_fourier_basis(
                omega_sidereal * t_k, self.fourier_sidereal_orders[param_str]
            )
            Amat = np.hstack((Amat, Amat_ps))
            keys.extend(self.get_fourier_sidereal_coeff_keys(param_str))

        if self.fit_velocity:
            # add the velocity terms
            Amat_dot = get_poly_basis_dt(
                t_k, self.poly_orders[param_str], dict["t_fit"], self.poly_types[param_str]
            )
            if self.fourier_orders[param_str] > 0:
                Amat_dot_p = get_fourier_basis_dt(
                    Phi, self.compute_Phidot(dict, t_k), self.fourier_orders[param_str]
                )
                Amat_dot = np.hstack((Amat_dot, Amat_dot_p))
            Amat = np.vstack((Amat, self.vel_scale * Amat_dot))
            residuals = np.hstack((residuals, self.vel_scale * residuals_dt))

        # Solve Optimization problem
        n = Amat.shape[0]
        x = cp.Variable(Amat.shape[1])

        if self.fit_obj == "lsq":
            # solve using scipy least squares
            coeffs, res_opt, rank, s = lstsq(Amat, residuals)

            res_rms = np.sqrt(res_opt) / n
            if param_str == "M":
                res_rms = np.mod(res_rms, np.pi * 2)
            if print_result:
                print("  [Fit {}] ".format(param_str))
                print("    Residual after fitting: ", res_rms)
                print("    Coefficients: ", coeffs)

        else:
            if self.fit_obj == "lsq-cvx":
                obj = cp.Minimize(cp.sum_squares(Amat @ x - residuals))
            elif self.fit_obj == "max":
                obj = cp.Minimize(cp.norm(Amat @ x - residuals, "inf"))
            elif self.fit_obj == "cvar":
                eta = cp.Variable()
                absres = cp.abs(Amat @ x - residuals)
                cvar = eta + (1 / (1 - self.cvar_beta)) * cp.sum(cp.pos(absres - eta)) / n
                obj = cp.Minimize(cvar)
            else:
                print("Invalid fit objective: ", self.fit_obj)
                return None

            const = []
            if fix_endpoints:
                const.append(Amat[0, :] @ x == residuals[0])
                const.append(Amat[-1, :] @ x == residuals[-1])
            prob = cp.Problem(obj, constraints=const)
            prob.solve(solver=cp.ECOS, max_iters=100000)
            coeffs = x.value

            if print_result:
                print("  [Fit {}] ".format(param_str))
                print("    Status                : ", prob.status)
                print("    Residual after fitting: ", prob.value)
                print("    Coefficients: ", coeffs)

        # replace the coefficients in the ephemeris
        ephem = self.replace_ephem(ephem, keys, coeffs)

        # plot the residuals and fit
        plot_fit = False

        if plot_fit:
            fitvals = Amat @ coeffs
            plt.figure(figsize=(10, 4))
            if param_str == "a":
                residuals = (residuals - 1.131e7) / 1000
                fitvals = (fitvals - 1.131e7) / 1000
            if param_str == "M":
                residuals = wrapToPi(residuals)
                fitvals = wrapToPi(fitvals)
            plt.plot(t_k + dict["t_ref"], residuals, "k.", label="data", markersize=4)
            plt.plot(
                t_k + dict["t_ref"],
                fitvals,
                "r-",
                label="fit",
                linewidth=2,
            )
            plt.xlabel("Time (s)")
            plt.ylabel("{} residual".format(param_str))
            plt.title("Fit of {}".format(param_str))
            plt.legend()
            plt.grid()
            plt.show()

        return ephem

    def fit_ruil_coeffs(self, t_data, param_data, dt_param_data, ephem, print_result=False):

        dict = self.ephem2dict(ephem)

        # extract the fixed parameters
        t_ref = dict["t_ref"]
        t_k = t_data - t_ref
        ndata = t_k.size

        a_0 = dict["a_0"]
        T_orbit = 2 * np.pi * np.sqrt(a_0**3 / self.GM)
        Phi_tmp = 2 * np.pi / T_orbit * t_k  # tempoaral variable for fourier

        # 1. fit a ----------------------------------------------------------------
        ephem = self.fit_coeffs(
            param_data["a"],
            "a",
            ephem,
            t_k,
            Phi_tmp,
            residuals_dt=dt_param_data["a"],
            print_result=print_result,
        )
        a_k = self.compute_param(self.ephem2dict(ephem), "a", t_k, Phi_tmp, 0)

        # 2. fit M ----------------------------------------------------------------
        n = np.sqrt(
            self.GM * np.ones_like(a_k) / a_k**3
        )  # Todo: replace this with integration (trapz?)

        res_M = np.unwrap(param_data["M"] - n * t_k)
        ephem = self.fit_coeffs(
            res_M,
            "M",
            ephem,
            t_k,
            Phi_tmp,
            print_result=print_result,
            residuals_dt=dt_param_data["M"] - n,
        )
        M_k = self.compute_param(self.ephem2dict(ephem), "M", t_k, Phi_tmp, n * t_k)  # mean anomaly

        # 3. fit e ----------------------------------------------------------------
        ephem = self.fit_coeffs(
            param_data["e"],
            "e",
            ephem,
            t_k,
            Phi_tmp,
            print_result=print_result,
            residuals_dt=dt_param_data["e"],
        )
        e_k = self.compute_param(self.ephem2dict(ephem), "e", t_k, Phi_tmp, 0)  # eccentricity

        # 4. solve Kepler's equation -----------------------------------------------------
        E_k = M_k
        for _ in range(5):
            E_k = E_k - (E_k - e_k * np.sin(E_k) - M_k) / (1 - e_k * np.cos(E_k))

        nu_k = 2 * np.arctan(np.sqrt((1 + e_k) / (1 - e_k)) * np.tan(E_k / 2))  # true anomaly

        # fit nu to improve convergence -------------------------------------------
        res_nu = np.unwrap(param_data["f"] - nu_k)
        ephem = self.fit_coeffs(
            res_nu,
            "f",
            ephem,
            t_k,
            Phi_tmp,
            print_result=print_result,
            residuals_dt=dt_param_data["f"],
        )
        nu_k = self.compute_param(self.ephem2dict(ephem), "f", t_k, Phi_tmp, 0)  # true anomaly

        # 5. compute argument of periapsis ------------------------------------------------
        ephem = self.fit_coeffs(
            param_data["w"],
            "w",
            ephem,
            t_k,
            Phi_tmp,
            print_result=print_result,
            residuals_dt=dt_param_data["w"],
        )
        w_k = self.compute_param(
            self.ephem2dict(ephem), "w", t_k, Phi_tmp, 0
        )  # argument of perigee

        # 5. compute the argument of latitude --------------------------------
        Phi = w_k + nu_k  # constant term

        # 6. fit radius -------------------------------------------------------------
        r0 = a_k * (1 - e_k * np.cos(E_k))
        r0_dot = self.compute_rdot(self.ephem2dict(ephem), t_k)
        ephem = self.fit_coeffs(
            param_data["r"] - r0,
            "r",
            ephem,
            t_k,
            Phi,
            print_result=print_result,
            fix_endpoints=False,
            residuals_dt=dt_param_data["r"] - r0_dot,
        )

        # 7. argument of latitude ---------------------------------------------
        res_u = np.unwrap(param_data["u"] - Phi)
        Phi_dot = self.compute_Phidot(self.ephem2dict(ephem), t_k)
        ephem = self.fit_coeffs(
            res_u,
            "u",
            ephem,
            t_k,
            Phi,
            print_result=print_result,
            fix_endpoints=False,
            residuals_dt=dt_param_data["u"] - Phi_dot,
        )

        # 3. inclination -------------------------------------------------
        ephem = self.fit_coeffs(
            param_data["i"],
            "i",
            ephem,
            t_k,
            Phi,
            print_result=print_result,
            fix_endpoints=False,
            residuals_dt=dt_param_data["i"],
        )

        # 4. longitude of ascending node ----------------------------------------
        res_lambda = np.unwrap(param_data["l"] + self.omega_b * t_k)
        ephem = self.fit_coeffs(
            res_lambda,
            "l",
            ephem,
            t_k,
            Phi,
            print_result=print_result,
            fix_endpoints=False,
            residuals_dt=dt_param_data["l"] + self.omega_b,
        )

        return ephem

    def init_guess(self, t_data, rvbf_data, print_result=False):
        """
        Initialize the coefficients for the keplarian ephemeris

        Args:
            t_data: array of times
            rv_data: array of position and velocity vectors in body frame
        """
        lent = t_data.size
        idx = int(lent / 2)  # use mid index
        t_ref = t_data[idx]
        rv_ref = rvbf_data[idx]

        t_fit = t_data[-1] - t_data[0]  # time span

        if print_result:
            print("  [Initial guess] ")
            print("    t_ref: ", t_ref)
            print("    t_fit: ", t_fit)

        coe = pnt.cart_to_classical(rv_ref, self.GM)
        dict = self.coe2dict(coe, t_ref, t_fit)

        if print_result:
            print(dict)

        # initialize the coefficients
        ephem = np.array(self.dict2ephem(dict))

        return t_ref, t_fit, ephem

    def generate_param_data(self, t_data, rv_data):

        # convert to ruil
        param_data = {
            "a": np.zeros(t_data.size),
            "e": np.zeros(t_data.size),
            "M": np.zeros(t_data.size),
            "f": np.zeros(t_data.size),
            "w": np.zeros(t_data.size),
            "r": np.zeros(t_data.size),
            "u": np.zeros(t_data.size),
            "i": np.zeros(t_data.size),
            "l": np.zeros(t_data.size),
        }
        lambda_prev = 0
        u_prev = 0
        M_prev = 0
        for i, t_k in enumerate(t_data):
            r = np.linalg.norm(rv_data[i, :3])
            coe = pnt.cart_to_classical(rv_data[i, :], self.GM)
            a = coe[0]
            e = coe[1]
            inc = coe[2]
            lamda = coe[3]
            w = coe[4]
            f = pnt.mean_to_true_anomaly(coe[5], coe[1])
            u = wrapTo2Pi(w + f)
            M = wrapTo2Pi(coe[5])

            # make sure the angles are continuous
            if i > 0:
                if abs(u - u_prev) > np.pi:
                    u += 2 * np.pi * np.sign(u_prev - u)
                if abs(lamda - lambda_prev) > np.pi:
                    lamda += 2 * np.pi * np.sign(lambda_prev - lamda)
                if abs(M - M_prev) > np.pi:
                    M += 2 * np.pi * np.sign(M_prev - M)

            u_prev = copy.deepcopy(u)
            lambda_prev = copy.deepcopy(lamda)
            M_prev = copy.deepcopy(M)

            param_data["a"][i] = a
            param_data["e"][i] = e
            param_data["M"][i] = M
            param_data["f"][i] = f
            param_data["w"][i] = w
            param_data["r"][i] = r
            param_data["u"][i] = u
            param_data["i"][i] = inc
            param_data["l"][i] = lamda

        dt_param_data = {}
        for key in param_data.keys():
            dt_param_data[key] = np.gradient(param_data[key], t_data)

        return param_data, dt_param_data

    def fit(
        self,
        t_data,
        rv_data,
        fit_obj="lsq",
        fit_velocity=False,
        vel_scale=1000,
        cvar_beta=0.95,
        print_result=False,
        use_scipy_lsq=True,
    ):
        """
        Fit the keplarian ephemeris to the data

        Args:
            t_data: array of times
            rv_data: array of position and velocity vectors in intertial frame attached to the body
        """
        self.fit_velocity = fit_velocity
        self.vel_scale = vel_scale
        self.fit_obj = fit_obj
        self.cvar_beta = cvar_beta
        self.use_scipy_lsq = use_scipy_lsq

        t_ref, t_fit, ephem = self.init_guess(t_data, rv_data, print_result=print_result)

        param_data, dt_param_data = self.generate_param_data(t_data, rv_data)

        # Optimize the coefficients for the angles and r seperately
        ephem = self.fit_ruil_coeffs(
            t_data, param_data, dt_param_data, ephem, print_result=print_result
        )

        return ephem

    def plot_fit_ruil(
        self,
        t_data,
        rvf_data,
        rvfw_data,
        ephem,
        plot_diff_pos=False,
        plot_diff_vel=False,
        plot_init=False,
        axes_ruil=None,
        axes_comp=None,
        plot_compare=False,
        moving_average_size=60,
    ):

        t_ref, t_fit, x0 = self.init_guess(t_data, rvf_data)
        ms = 1  # marker size
        lw = 1  # line width

        lent = rvf_data.shape[0]
        ruil0 = np.zeros((lent, 4))
        ruil = np.zeros((lent, 4))
        ruil_true = np.zeros((lent, 4))
        dt_ruil0 = np.zeros((lent, 4))
        dt_ruil = np.zeros((lent, 4))
        dt_ruil_true = np.zeros((lent, 4))

        # compute fitted trajectory and errors --------------------------------
        ruil0, dt_ruil0 = self.ephem2ruil(t_data, x0, compute_velocity=True)
        ruil, dt_ruil = self.ephem2ruil(t_data, ephem, compute_velocity=True)

        # convert to numpy array ---------------------------------------------
        ruil0 = np.array(ruil0)
        ruil = np.array(ruil)
        dt_ruil0 = np.array(dt_ruil0)
        dt_ruil = np.array(dt_ruil)

        r = np.linalg.norm(rvf_data[:, :3], axis=1)
        coe = pnt.cart_to_classical(rvf_data, self.GM)
        a = coe[:, 0]
        e = coe[:, 1]
        p = a * (1 - e**2)  # semi-latus rectum
        w = coe[:, 4]
        inc = coe[:, 2]
        lamda = coe[:, 3]
        f = pnt.mean_to_true_anomaly(coe[:, 5], coe[:, 1])
        u = w + f
        ruil_true = np.vstack([r, u, inc, lamda]).T  # T x 4

        # fix sudden jumps in u and l --------------------------------
        idxs = [1, 3]
        for idx in idxs:
            angle = ruil[:, idx]
            angle_true = ruil_true[:, idx]
            angle_init = ruil0[:, idx]

            for i in range(1, lent):
                if angle[i] - angle[i - 1] > np.pi:
                    ruil[i:, idx] -= 2 * np.pi
                elif angle[i] - angle[i - 1] < -np.pi:
                    ruil[i:, idx] += 2 * np.pi

                if angle_true[i] - angle_true[i - 1] > np.pi:
                    ruil_true[i:, idx] -= 2 * np.pi
                elif angle_true[i] - angle_true[i - 1] < -np.pi:
                    ruil_true[i:, idx] += 2 * np.pi

                if angle_init[i] - angle_init[i - 1] > np.pi:
                    ruil0[i:, idx] -= 2 * np.pi
                elif angle_init[i] - angle_init[i - 1] < -np.pi:
                    ruil0[i:, idx] += 2 * np.pi

        # for true velocity terms, compute by differentiation ------------------
        for i in range(4):
            if i == 0:
                dt_ruil_true[:, i] = np.sqrt(self.GM * np.ones_like(p) / p) * e * np.sin(f)
            else:
                # dt_ruil_true[:-1, i] = np.diff(ruil_true[:, i])/np.diff(t_data)
                dt_ruil_true[:, i] = np.gradient(ruil_true[:, i], t_data)

        t_data = (t_data - t_data[0]) / 60

        # conversion of units --------------------------------
        # convert to degrees
        for i in range(3):
            ruil0[:, i + 1] = np.rad2deg(ruil0[:, i + 1])
            ruil[:, i + 1] = np.rad2deg(ruil[:, i + 1])
            ruil_true[:, i + 1] = np.rad2deg(ruil_true[:, i + 1])

            dt_ruil0[:, i + 1] = np.rad2deg(dt_ruil0[:, i + 1])
            dt_ruil[:, i + 1] = np.rad2deg(dt_ruil[:, i + 1])

        # scale r to meters
        # ruil0[:, 0] = ruil0[:, 0] * 1e3
        # ruil[:, 0] = ruil[:, 0] * 1e3
        # ruil_true[:, 0] = ruil_true[:, 0] * 1e3
        # dt_ruil0[:, 0] = dt_ruil0[:, 0] * 1e3
        # dt_ruil[:, 0] = dt_ruil[:, 0] * 1e3
        # dt_ruil_true[:, 0] = dt_ruil_true[:, 0] * 1e3

        # plot r, u, i ---------------------------------------------------------
        if axes_ruil is None:
            fig, axes = plt.subplots(4, 2, figsize=(12, 12))
        else:
            axes = axes_ruil

        if plot_diff_pos:
            labels = [
                r"$\Delta r$ [m]",
                r"$\Delta u$ [deg]",
                r"$\Delta i$ [deg]",
                r"$\Delta \lambda$ [deg]",
            ]
        else:
            labels = ["$r$ [k]", "$u$ [deg]", "$i$ [deg]", r"$\lambda$ [deg]"]

        masize = moving_average_size
        for i, ax in enumerate(axes[:, 0]):
            if plot_diff_pos:
                if plot_init:
                    ax.plot(
                        t_data,
                        (ruil0[:, i] - ruil_true[:, i]),
                        "ro--",
                        markersize=ms,
                        linewidth=lw,
                        label="True - Initial",
                    )
                ax.plot(
                    t_data,
                    (ruil[:, i] - ruil_true[:, i]),
                    "bo-",
                    markersize=ms,
                    linewidth=lw,
                    label="Fit - True",
                )
                # plot moving average
                maplot = np.convolve(
                    ruil[:, i] - ruil_true[:, i], np.ones(masize) / masize, mode="same"
                )
                ax.plot(
                    t_data[masize : lent - masize],
                    maplot[masize : lent - masize],
                    "m-",
                    markersize=ms,
                    linewidth=lw * 2,
                    label="Fit - True (MA)",
                )
            else:
                ax.plot(t_data, ruil_true[:, i], "ko-", markersize=ms, linewidth=lw, label="True")
                if plot_init:
                    ax.plot(
                        t_data, ruil0[:, i], "ro--", markersize=ms, linewidth=lw, label="Initial"
                    )
                ax.plot(t_data, ruil[:, i], "bo-", markersize=ms, linewidth=lw, label="Fit")
            ax.set_xlabel("Time [min]")
            ax.set_ylabel(labels[i])
            ax.legend()
            # if i > 0 and plot_diff_pos:
            #     # ax.set_ylim(-1e-4, 1e-4)
            #     # ax.set_yticks(np.arange(-10e-5, 11e-5, 2.5e-5))
            #     # ax.set_yticklabels(['-10.0e-5', '-7.5e-5', '-5.0e-5', '-2.5e-5', '0', '2.5e-5', '5.0e-5', '7.5e-5', '10.0e-5'])
            ax.grid(True)

        # dr, du, di, dlambda
        if plot_diff_vel:
            labels = [
                r"$\Delta \dot{r}$ [m/s]",
                r"$\Delta \dot{u}$ [deg/s]",
                r"$\Delta \dot{i}$ [deg/s]",
                r"$\Delta \dot{\lambda}$ [deg/s]",
            ]
        else:
            labels = [
                r"$\dot{r}$ [m/s]",
                r"$\dot{u}$ [deg/s]",
                r"$\dot{i}$ [deg/s]",
                r"$\dot{\lambda}$ [deg/s]",
            ]

        # remove first and last point
        t_data = t_data[1:-1]
        dt_ruil0 = dt_ruil0[1:-1]
        dt_ruil = dt_ruil[1:-1]
        dt_ruil_true = dt_ruil_true[1:-1]

        for i, ax in enumerate(axes[:, 1]):
            if plot_diff_vel:
                if plot_init:
                    ax.plot(
                        t_data,
                        (dt_ruil0[:, i] - dt_ruil_true[:, i]),
                        "ro--",
                        markersize=ms,
                        linewidth=lw,
                        label="True - Initial",
                    )
                ax.plot(
                    t_data,
                    (dt_ruil[:, i] - dt_ruil_true[:, i]),
                    "bo-",
                    markersize=ms,
                    linewidth=lw,
                    label="Fit - True",
                )
                # plot moving average
                maplot = np.convolve(
                    dt_ruil[:, i] - dt_ruil_true[:, i], np.ones(masize) / masize, mode="same"
                )
                ax.plot(
                    t_data[masize : lent - masize],
                    maplot[masize : lent - masize],
                    "m-",
                    markersize=ms,
                    linewidth=lw * 2,
                    label="Fit - True (MA)",
                )
            else:
                ax.plot(
                    t_data, dt_ruil_true[:, i], "ko-", markersize=ms, linewidth=lw, label="True"
                )
                if plot_init:
                    ax.plot(
                        t_data, dt_ruil0[:, i], "ro--", markersize=ms, linewidth=lw, label="Initial"
                    )
                ax.plot(t_data, dt_ruil[:, i], "b-", markersize=ms, linewidth=lw, label="Fit")
            ax.set_xlabel("Time [min]")
            ax.set_ylabel(labels[i])
            ax.legend()
            # ax.set_xlim([3, 4])
            ax.grid(True)

        plt.tight_layout()

        return axes
