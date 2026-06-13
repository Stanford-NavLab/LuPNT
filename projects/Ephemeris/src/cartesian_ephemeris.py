import numpy as np
import pylupnt as pnt
import cvxpy as cp
import matplotlib.pyplot as plt
from copy import deepcopy

try:
    from .ephemeris import Ephemeris
    from .orbit_utils import wrapToPi, wrapTo2Pi, rms
    from .orbit_utils import cart_to_mqoe, mqoe_to_cart, mqoe_to_coe, coe_to_mqoe
    from .basis import get_poly_basis, get_poly_basis_dt, get_fourier_basis, get_fourier_basis_dt
except ImportError:
    from ephemeris import Ephemeris
    from orbit_utils import wrapToPi, wrapTo2Pi, rms
    from orbit_utils import cart_to_mqoe, mqoe_to_cart, mqoe_to_coe, coe_to_mqoe
    from basis import get_poly_basis, get_poly_basis_dt, get_fourier_basis, get_fourier_basis_dt


class CartesianEphemeris(Ephemeris):
    def __init__(
        self,
        order,
        use_kep=True,
        use_rsw=False,
        use_fourier=False,
        use_meq=False,
        poly_type="chebyshev",
        convert_to_coe=False,
        body=pnt.MOON,
        print_info=False,
    ):
        super().__init__(body)
        self.order = order
        self.use_kep = use_kep
        self.use_rsw = use_rsw
        self.use_fourier = use_fourier
        self.use_meq = use_meq
        self.poly_type = poly_type
        self.convert_to_coe = convert_to_coe

        if not self.use_kep and self.use_rsw:
            print(
                "Warning: Cannot use RSW coordinates without Keplerian elements. Set use_kep=True."
            )
            self.use_kep = True

        if not self.use_kep and self.use_fourier:
            print(
                "Warning: Cannot use Fourier coefficients without Keplerian elements. Set use_kep=True."
            )
            self.use_kep = True

        self.create_ephem_dict(print_info)

    def create_ephem_dict(self, print_info=False):

        self.idx_keys = {}
        self.idx_keys["t_ref"] = 0
        self.idx_keys["t_fit"] = 1

        self.keys_list = ["t_ref", "t_fit"]

        if self.use_kep:
            if self.use_meq:  # modified equinoctial elements
                self.oe_params = ["p_0", "f_0", "g_0", "h_0", "k_0", "L_0"]
            else:
                self.oe_params = ["A_0", "e_0", "i_0", "l_0", "w_0", "M_0"]
            for ephem in self.oe_params:
                self.idx_keys[ephem] = len(self.idx_keys)

            self.keys_list.extend(self.oe_params)

        if self.use_rsw:
            self.poly_coeff_str = ["a", "c", "r"]
        else:
            self.poly_coeff_str = ["x", "y", "z"]

        for i, cstr in enumerate(self.poly_coeff_str):
            for j in range(self.order + 1):
                key = f"{cstr}_{j}"
                self.idx_keys[key] = len(self.idx_keys)
                self.keys_list.append(key)

        if self.use_fourier:
            for i, cstr in enumerate(self.poly_coeff_str):
                key = f"C_{cstr}_{self.order}"  # cosine
                self.idx_keys[key] = len(self.idx_keys)
                self.keys_list.append(key)

                key = f"S_{cstr}_{self.order}"  # sine
                self.idx_keys[key] = len(self.idx_keys)
                self.keys_list.append(key)

        self.n_params = len(self.keys_list)

        if print_info:
            print("key list: ", self.keys_list)
            print("Use Keplerian elements: ", self.use_kep)
            print("Use RSW coordinates: ", self.use_rsw)
            print("Chebyshev Order: ", self.order)
            print(f"Number of parameters: {self.n_params}")
        # for key, value in self.idx_keys.items():
        #     print(f"{key}: {value}")

    def oe2dict(self, coe, t_ref, t_fit):
        dict = {}

        for i, key in enumerate(self.keys_list):
            if key == "t_ref":
                dict[key] = t_ref
            elif key == "t_fit":
                dict[key] = t_fit
            elif key in self.oe_params:
                j = self.oe_params.index(key)
                dict[key] = coe[j]
            else:
                dict[key] = 0

        return dict

    def print_ephem(self, ephem):
        ephem_dict = self.ephem2dict(ephem)
        coeffs = np.zeros((3, self.order + 1))
        for key, value in ephem_dict.items():
            split_key = key.split("_")
            if split_key[0] in self.poly_coeff_str:
                if split_key[0] == self.poly_coeff_str[0]:
                    idx = 0
                elif split_key[0] == self.poly_coeff_str[1]:
                    idx = 1
                elif split_key[0] == self.poly_coeff_str[2]:
                    idx = 2
                else:
                    print("Key not found in poly_coeff_str: ", split_key[0])
                coeffs[idx, int(split_key[1])] = value
            else:
                if key in ["i_0", "w_0", "l_0", "M_0", "L_0"]:
                    value = np.rad2deg(value)
                if abs(value) >= 0.01:
                    print("{0:}: {1:.4f}".format(key, value))
                else:
                    print("{0:}: {1:.2e}".format(key, value))

        for i, cstr in enumerate(self.poly_coeff_str):
            print(f"{cstr} coeffs: ", coeffs[i, :])

    def coe_ephem2cartu(self, t_k, dict, compute_velocity=False):

        a = dict["A_0"]
        e = dict["e_0"]
        i_k = dict["i_0"]
        l_0 = dict["l_0"]
        w = dict["w_0"]
        M0 = dict["M_0"]

        # if t_k is a scalar, convert it to a 1D array
        coe = np.array([a, e, i_k, l_0, w, M0])

        return self.coe_vec2cartu(t_k, coe, compute_velocity)

    def coe_vec2cartu(self, t_k, coe, compute_velocity=False):

        a, e, i_k, l_0, w, M0 = coe

        n = np.sqrt(self.GM * np.ones_like(a) / a**3)

        # Compute the mean anomaly
        M_k = M0 + n * t_k

        # Solve Kepler's equation for E
        E = M_k
        for _ in range(5):
            E = E - (E - e * np.sin(E) - M_k) / (1 - e * np.cos(E))

        # Compute the true anomaly
        nu_k = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))

        # Compute the semi-latus rectum
        u_k = w + nu_k

        # Compute the radius
        r_k = a * (1 - e * np.cos(E))

        lambda_k = l_0 - self.omega_b * t_k

        x_k = r_k * np.cos(u_k)
        y_k = r_k * np.sin(u_k)
        x_pos = x_k * np.cos(lambda_k) - y_k * np.cos(i_k) * np.sin(lambda_k)
        y_pos = x_k * np.sin(lambda_k) + y_k * np.cos(i_k) * np.cos(lambda_k)
        z_pos = y_k * np.sin(i_k)

        if compute_velocity:
            Edot = n / (1 - e * np.cos(E))
            udot = Edot * np.sqrt(1 - e**2) / (1 - e * np.cos(E))
            rdot = a * e * np.sin(E) * Edot
            lambdadot = -self.omega_b
            didt = 0
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

            params = {
                "r_k": r_k,
                "u_k": u_k,
                "lambda_k": lambda_k,
                "i_k": i_k * np.ones_like(u_k),
                "M_k": M_k,
                "rdot": rdot,
                "udot": udot,
                "lambdadot": lambdadot,
                "didt": didt,
            }

            return np.vstack([x_pos, y_pos, z_pos, x_vel, y_vel, z_vel]).T, u_k, udot, params

        else:
            params = {"r_k": r_k, "u_k": u_k, "lambda_k": lambda_k, "i_k": i_k, "M_k": M_k}
            return np.vstack([x_pos, y_pos, z_pos]).T, u_k, None, params

    def mqoe_ephem2cartu(self, t_k, dict, compute_velocity=False):
        p0 = dict["p_0"]
        f0 = dict["f_0"]
        g0 = dict["g_0"]
        h0 = dict["h_0"]
        k0 = dict["k_0"]
        L0 = dict["L_0"]

        qoe = np.array([p0, f0, g0, h0, k0, L0])

        if self.convert_to_coe:
            # first convert to coe
            coe = mqoe_to_coe(qoe, use_true_anom=False)
            res, u_k, udot_k, params_coe = self.coe_vec2cartu(
                t_k, coe, compute_velocity=compute_velocity
            )

            if t_k.ndim == 0:
                lam_k = params_coe["lambda_k"]
                M_k = params_coe["M_k"]
                coe_vec = np.array([coe[0], coe[1], coe[2], lam_k, coe[4], M_k])
            else:
                coe_vec = np.zeros((t_k.size, 6))
                coe_vec[:, 0] = coe[0]  # a
                coe_vec[:, 1] = coe[1]  # e
                coe_vec[:, 2] = coe[2]  # i
                coe_vec[:, 3] = params_coe["lambda_k"]  # l
                coe_vec[:, 4] = coe[4]  # w
                coe_vec[:, 5] = params_coe["M_k"]  # M

            qoe_vec = coe_to_mqoe(coe_vec)

            params = {}
            if t_k.ndim == 0:
                params["p_k"] = qoe_vec[0]
                params["f_k"] = qoe_vec[1]
                params["g_k"] = qoe_vec[2]
                params["h_k"] = qoe_vec[3]
                params["k_k"] = qoe_vec[4]
                params["L_k"] = qoe_vec[5]
            else:
                params["p_k"] = qoe_vec[:, 0]
                params["f_k"] = qoe_vec[:, 1]
                params["g_k"] = qoe_vec[:, 2]
                params["h_k"] = qoe_vec[:, 3]
                params["k_k"] = qoe_vec[:, 4]
                params["L_k"] = qoe_vec[:, 5]

            return res, u_k, udot_k, params

        else:
            return self.mqoe_vec2cartu(t_k, qoe, compute_velocity)

    def mqoe_vec2cartu(self, t_k, qoe, compute_velocity=False):

        p0, f0, g0, h0, k0, L0 = qoe

        e0 = np.sqrt(f0**2 + g0**2)
        a0 = p0 / (1 - e0**2)

        # correct for Omega
        Omega0 = np.arctan2(k0, h0)
        Omega_dot = -self.omega_b
        Omega_k = Omega0 + Omega_dot * t_k
        Omega_w0 = np.arctan2(g0, f0)  # Omega + w
        Omega_w_k = Omega_w0 + Omega_dot * t_k  # Omega + w

        f_k = e0 * np.cos(Omega_w_k)
        g_k = e0 * np.sin(Omega_w_k)
        tan_i20 = np.sqrt(h0**2 + k0**2)
        h_k = tan_i20 * np.cos(Omega_k)
        k_k = tan_i20 * np.sin(Omega_k)

        # compute angles
        nu0 = L0 - Omega_w0
        E0 = 2 * np.arctan(np.sqrt((1 - e0) / (1 + e0)) * np.tan(nu0 / 2))
        M0 = E0 - e0 * np.sin(E0)
        n = np.sqrt(self.GM / a0**3)
        M_k = n * t_k + M0
        # Solve Kepler's equation for E
        E_k = M_k
        for _ in range(5):
            E_k = E_k - (E_k - e0 * np.sin(E_k) - M_k) / (1 - e0 * np.cos(E_k))
        # Convert to nu
        nu_k = 2 * np.arctan(np.sqrt((1 + e0) / (1 - e0)) * np.tan(E_k / 2))

        # Compute L_k
        L_k = wrapTo2Pi(L0 + (nu_k - nu0) + Omega_dot * t_k)
        cLk = np.cos(L_k)
        sLk = np.sin(L_k)

        w_k = 1 + f_k * np.cos(L_k) + g_k * np.sin(L_k)  # 1 + e * cos(nu)
        r = p0 / w_k
        s2 = 1 + h_k**2 + k_k**2
        alpha2 = h_k**2 - k_k**2

        r_s2 = r / s2
        thk = 2 * k_k * h_k

        hsl_kcl = h_k * sLk - k_k * cLk
        hcl_ksl = h_k * cLk + k_k * sLk

        u_k = np.arctan2(hsl_kcl, hcl_ksl)
        Edot = n / (1 - e0 * np.cos(E_k))
        udot_k = Edot * np.sqrt(1 - e0**2) / (1 - e0 * np.cos(E_k))

        x_k_tmp = (1 + alpha2) * cLk + thk * sLk
        y_k_tmp = (1 - alpha2) * sLk + thk * cLk
        x_k = r_s2 * x_k_tmp
        y_k = r_s2 * y_k_tmp
        z_k = 2 * r_s2 * (hsl_kcl)

        params = {
            "p_k": p0 * np.ones_like(f_k),
            "f_k": f_k,
            "g_k": g_k,
            "h_k": h_k,
            "k_k": k_k,
            "L_k": L_k,
        }

        if compute_velocity:
            # v = np.sqrt(self.GM/p0) / s2
            nudot = np.sqrt(self.GM * p0) * (w_k / p0) ** 2
            Ldot = nudot + Omega_dot
            fdot = -e0 * np.sin(Omega_w_k) * Omega_dot
            gdot = e0 * np.cos(Omega_w_k) * Omega_dot
            wdot = (
                fdot * np.cos(L_k)
                + gdot * np.sin(L_k)
                - f_k * Ldot * np.sin(L_k)
                + g_k * Ldot * np.cos(L_k)
            )
            rdot_s2 = (p0 / w_k**2) * wdot / s2
            alpha2_dot = (
                -(tan_i20**2) * 2 * np.sin(2 * Omega_k) * Omega_dot
            )  # alpha2 = tan(i/2)^2 * cos(2*Omega)
            thk_dot = (
                (tan_i20) * 2 * np.cos(2 * Omega_k) * Omega_dot
            )  # 2hk = tan(i/2) * sin(2*Omega)
            hk_dot = -tan_i20 * np.sin(Omega_k) * Omega_dot
            kk_dot = tan_i20 * np.cos(Omega_k) * Omega_dot

            x_k_tmp_dot = Ldot * (
                -(1 + alpha2) * sLk + alpha2_dot * cLk + thk * cLk + thk_dot * sLk
            )
            y_k_tmp_dot = Ldot * ((1 - alpha2) * cLk - alpha2_dot * sLk - thk * sLk + thk_dot * cLk)
            hsl_kcl_dot = Ldot * (h_k * cLk + hk_dot * sLk - (kk_dot * cLk - k_k * sLk))

            vx_k = rdot_s2 * x_k_tmp + r_s2 * x_k_tmp_dot
            vy_k = rdot_s2 * y_k_tmp + r_s2 * y_k_tmp_dot
            vz_k = 2 * rdot_s2 * hsl_kcl + 2 * r_s2 * hsl_kcl_dot

            # vx_k = rdot_s2 * (-sLk - alpha2 * sLk + thk * cLk - g_k + f_k * thk - alpha2 * g_k)
            # vy_k = rdot_s2 * (cLk - alpha2 * cLk - thk * sLk + f_k - g_k * thk - alpha2 * f_k)
            # vz_k = 2*rdot_s2 * (h_k * cLk + k_k * sLk + f_k * h_k + g_k * k_k)

            return np.vstack([x_k, y_k, z_k, vx_k, vy_k, vz_k]).T, u_k, udot_k, params

        else:
            return np.vstack([x_k, y_k, z_k]).T, u_k, None, params

    def poly2cart(self, t_k, dict, compute_velocity=False, u_k=None, udot_k=None):
        """
        Convert Chebyshev coefficients to Cartesian coordinates.

        Parameters
        ----------
        t_k : array_like
            Time in seconds since the reference time.
        dict : dict
            Dictionary containing Chebyshev coefficients.
        compute_velocity : bool, optional
            If True, compute velocity as well. Default is False.

        Returns
        -------
        res : ndarray
            Cartesian coordinates (and velocity if compute_velocity is True).
        """
        # normalized time
        # t_k = -t_fit/2, ..., 0, ..., t_fit/2
        # to
        # t_k = -1, ..., 0, ..., 1
        # t_k_center = (t_k[0] + t_k[-1]) / 2
        # t_k = t_k - t_k_center  # bring it to zero center
        z_k = 2 * t_k / dict["t_fit"]  # 1 x T
        lent = t_k.size

        # Chebyshev coefficients
        T = get_poly_basis(t_k, self.order, dict["t_fit"], self.poly_type)

        if self.use_fourier:
            Tf = get_fourier_basis(u_k, 1)
            T = np.hstack((T, Tf))

        coeffs = np.zeros((T.shape[1], 3))
        for i in range(self.order + 1):
            if self.use_rsw:
                coeffs[i, 0] = dict[f"a_{i}"]  # along-track
                coeffs[i, 1] = dict[f"c_{i}"]  # cross-track
                coeffs[i, 2] = dict[f"r_{i}"]  # radial
            else:
                coeffs[i, 0] = dict[f"x_{i}"]
                coeffs[i, 1] = dict[f"y_{i}"]
                coeffs[i, 2] = dict[f"z_{i}"]

        if self.use_fourier:
            for i, cstr in enumerate(self.poly_coeff_str):
                coeffs[self.order + 1, i] = dict[f"C_{cstr}_{self.order}"]
                coeffs[self.order + 2, i] = dict[f"S_{cstr}_{self.order}"]

        # (T x n_order) x (n_order x 3) = T x 3
        res_poly = T @ coeffs

        if compute_velocity:
            Tdot = get_poly_basis_dt(t_k, self.order, dict["t_fit"], self.poly_type)

            if self.use_fourier:
                Tdotf = get_fourier_basis_dt(u_k, udot_k, 1)
                Tdot = np.hstack((Tdot, Tdotf))

            res_poly_dot = Tdot @ coeffs  # T x 3
            res_poly = np.hstack((res_poly, res_poly_dot))  # T x 6

        return res_poly

    def rot_i2rsw(self, rv_coe):
        # if dimension is 1, add a dimension
        # inertial to RSW frame
        e_a = rv_coe[3:] / np.linalg.norm(rv_coe[3:])
        e_c = np.cross(rv_coe[:3], rv_coe[3:])
        e_c = e_c / np.linalg.norm(e_c)
        e_r = np.cross(e_a, e_c)
        e_r = e_r / np.linalg.norm(e_r)

        return np.vstack([e_a, e_c, e_r]).T  # 3 x 3

    def ephem2cart(self, t, ephem, compute_velocity=False, return_params=False, scale=None):

        dict = self.ephem2dict(ephem, scale=scale)
        t_ref = dict["t_ref"]

        t_k = t - t_ref
        lent = t.shape[0]
        params = {}

        if self.use_kep:
            if self.use_meq:
                rv_coe, u_k, udot_k, params = self.mqoe_ephem2cartu(
                    t_k, dict, compute_velocity=True
                )  # T x 6
            else:
                rv_coe, u_k, udot_k, params = self.coe_ephem2cartu(
                    t_k, dict, compute_velocity=True
                )  # T x 6
        else:
            rv_coe = np.zeros((lent, 6))
            u_k = None
            udot_k = None

        res_poly = self.poly2cart(t_k, dict, compute_velocity, u_k, udot_k)  # T x 3 or T x 6

        if not self.use_rsw:
            if compute_velocity:
                res = rv_coe + res_poly
            else:
                res = rv_coe[:, :3] + res_poly
        else:
            res = np.zeros((lent, 3 + int(3 * compute_velocity)))
            for ti in range(lent):
                M_i2rsw = self.rot_i2rsw(rv_coe[ti, :])  # 3 x 3
                M_rsw2i = M_i2rsw.T  # 3 x 3
                if compute_velocity:
                    Mrotrv = np.zeros((6, 6))
                    Mrotrv[:3, :3] = M_rsw2i
                    Mrotrv[3:, 3:] = M_rsw2i
                    res_poly_xyz = Mrotrv @ res_poly[ti]
                    res[ti] = rv_coe[ti] + res_poly_xyz  # T x 6
                else:
                    res_poly_xyz = M_rsw2i @ res_poly[ti]
                    res[ti] = rv_coe[ti, :3] + res_poly_xyz  # T x 3

        if return_params:
            return res, params
        else:
            return res

    def init_guess(self, t_data, rvbf_data):
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

        if self.use_kep:
            if self.use_meq:
                oe = cart_to_mqoe(rv_ref, self.GM)
                # coe = pnt.cart_to_classical(rv_ref, self.GM)
                # oe = coe_to_mqoe(coe)
            else:
                oe = pnt.cart_to_classical(rv_ref, self.GM)

            dict = self.oe2dict(oe, t_ref, t_fit)

            # initialize the coefficients
            ephem = self.dict2ephem(dict)
        else:
            ephem = np.zeros(self.n_params)
            ephem[0] = t_ref
            ephem[1] = t_fit

        return t_ref, t_fit, ephem

    def poly_fit(
        self,
        t_k,
        res_r,
        ephem,
        print_result=False,
        u_k=None,
        constrain_endpoints=False,
        L2_alpha=None,
    ):
        res_r = res_r * 1e-3  # convert to km

        # normalized time
        # t_k = -t_fit/2, ..., 0, ..., t_fit/2
        # to
        # t_k = -1, ..., 0, ..., 1
        dict = self.ephem2dict(ephem)

        # t_k_center = (t_k[0] + t_k[-1]) / 2
        # t_k = t_k - t_k_center  # bring it to zero center
        z_k = 2 * t_k / dict["t_fit"]  # 1 x T
        lent = t_k.size

        # Chebyshev coefficients
        T = get_poly_basis(t_k, self.order, dict["t_fit"], self.poly_type)

        # Fourier coefficients
        if self.use_fourier:
            Tf = get_fourier_basis(u_k, 1)
            T = np.hstack((T, Tf))

        res = [None] * 3
        if self.use_rsw:
            param_str = ["a", "c", "r"]
            labels = ["[along-track]", "[cross-track]", "[radial]"]
        else:
            param_str = ["x", "y", "z"]
            labels = ["[x]", "[y]", "[z]"]

        eps = 1e-9  # small value to avoid division by zero
        xtol = 1e-9  # tolerance for convergence
        n_iter_log = 30  # number of iterations for log regularization

        for i in range(3):
            # optimization problem
            # minimize ||Tx - res_r||_2
            n = T.shape[1]
            x = cp.Variable(n)

            if self.fit_obj == "lsq":
                y = res_r[:, i]
                x_sol, res_opt, rank, s = np.linalg.lstsq(T, y)
                res_rns = rms(res_opt)
                if print_result:
                    print(labels[i])
                    print("  Residual after fitting: ", res_opt)
                    print("  Polynomial coefficients: ", x_sol[: (self.order + 1)])
                    if self.use_fourier:
                        print("  Fourier coefficients   : ", x_sol[(self.order + 1) :])
                    start_res = T[0, :] @ x_sol - res_r[0, i]
                    end_res = T[-1, :] @ x_sol - res_r[-1, i]
                    print("  Start residual        : ", start_res)
                    print("  End residual          : ", end_res)

                # if L2_alpha is not None and L2_alpha > 0:
                #     # add L2 regularization
                #     for k in range(n_iter_log):
                #         weights = L2_alpha * np.ones_like(x_sol) / (np.abs(x_sol) + eps)
                #         T_new = np.vstack((T, np.sqrt(weights) * np.eye(n)))
                #         y_new = np.hstack((y, np.zeros(n)))
                #         x_sol, res_opt, rank, s = np.linalg.lstsq(T_new, y_new)

                # if use_l2:
                #     print("  ")
                #     print("  L2 regularization     : ", L2_alpha)
                #     print("  Residual norm         : ", res_rns)
                #     print("  L2 norm of coefficients: ", np.linalg.norm(x_sol, 2))

            else:
                if self.fit_obj == "lsq-cvx":
                    obj = cp.Minimize(cp.sum_squares(T @ x - res_r[:, i]))
                elif self.fit_obj == "max":
                    obj = cp.Minimize(cp.norm(T @ x - res_r[:, i], "inf"))
                elif self.fit_obj == "cvar":
                    eta = cp.Variable()
                    absres = cp.abs(T @ x - res_r[:, i])
                    cvar = eta + (1 / (1 - self.cvar_beta)) * cp.sum(cp.pos(absres - eta)) / lent
                    obj = cp.Minimize(cvar)
                elif self.fit_obj == "l1":
                    obj = cp.Minimize(cp.sum_squares(T @ x - res_r[:, i]) + L2_alpha * cp.norm1(x))
                elif self.fit_obj == "l1-log":
                    # first solve the l1 problem
                    obj = cp.Minimize(cp.sum_squares(T @ x - res_r[:, i]))
                else:
                    raise ValueError("fit_obj must be 'lsq' or 'max'")

                prob = cp.Problem(obj, [])
                prob.solve(solver=cp.ECOS, max_iters=1000)
                x_sol = x.value

                # iteratively refine the solution with log regularization
                if self.fit_obj == "l1-log":
                    for k in range(n_iter_log):
                        x = cp.Variable(n)
                        W = np.ones_like(x_sol) / (np.abs(x_sol) + eps)
                        W = np.diag(W)
                        obj = cp.Minimize(
                            cp.sum_squares(T @ x - res_r[:, i]) + L2_alpha * cp.sum(W @ cp.abs(x))
                        )
                        prob = cp.Problem(obj, [])
                        prob.solve(solver=cp.ECOS, max_iters=1000)
                        dx = np.linalg.norm(x.value - x_sol)
                        # print("  Iteration {}, change in x: {}".format(k + 1, dx))
                        if dx < xtol:
                            # print("Converged after {} iterations".format(k+1))
                            x_sol = x.value
                            break
                        x_sol = x.value

                if print_result:
                    print(labels[i])
                    print("  Status                : ", prob.status)
                    print("  Residual after fitting: ", prob.value)
                    print("  Polynomial coefficients: ", x_sol[: (n + 1)])
                    if self.use_fourier:
                        print("  Fourier coefficients   : ", x_sol[(n + 1) :])
                    start_res = T[0, :] @ x_sol - res_r[0, i]
                    end_res = T[-1, :] @ x_sol - res_r[-1, i]
                    print("  Start residual        : ", start_res)
                    print("  End residual          : ", end_res)

                if x_sol is None:
                    # resolve with least squares
                    y = res_r[:, i]
                    x_sol, res_opt, rank, s = np.linalg.lstsq(T, y)

            # store results (chebyshev)
            x_sol = x_sol * 1e3  # back to meters

            for k in range(self.order + 1):
                dict[f"{param_str[i]}_{k}"] = x_sol[k]

            if self.use_fourier:
                dict[f"C_{param_str[i]}_{k}"] = x_sol[self.order + 1]
                dict[f"S_{param_str[i]}_{k}"] = x_sol[self.order + 2]

            plot_result = False

            if plot_result:
                res[i] = T @ (x_sol * 1e-3) - res_r[:, i]
                plt.figure(figsize=(8, 3))
                plt.plot(t_k / 86400, res[i], "b-")
                plt.xlabel("Time [days]")
                plt.ylabel("Residual " + labels[i] + " [km]")
                plt.title("Residuals after fitting " + labels[i])
                plt.grid()
                plt.show()

        ephem = self.dict2ephem(dict)

        return ephem

    def fit(
        self, t_data, rv_data, fit_obj="lsq", print_result=False, cvar_beta=0.95, L2_alpha=None
    ):
        """
        Fit the ephemeris to the data

        Args:
            t_data: array of times
            rv_data: array of position and velocity vectors in body frame [N x 6]
            fit_obj: fitting objective, 'lsq' for least squares, 'max' for min-max, 'cvar' for conditional value at risk
            print_result: if True, print the fitting results
            cvar_beta: beta parameter for cvar fitting (only used if fit_obj='cvar')
            L2_alpha: L2 regularization parameter (only used if fit_obj='l1' or 'l1-log')
        Returns:
            ephem: fitted ephemeris coefficients

        """

        self.cvar_beta = cvar_beta

        t_ref, t_fit, ephem = self.init_guess(t_data, rv_data)
        t_k = t_data - t_ref
        self.fit_obj = fit_obj

        # compute the residuals
        if self.use_kep:
            if self.use_meq:
                rv_oe, u_k, _, _ = self.mqoe_ephem2cartu(
                    t_k, self.ephem2dict(ephem), compute_velocity=True
                )
            else:
                rv_oe, u_k, _, _ = self.coe_ephem2cartu(
                    t_k, self.ephem2dict(ephem), compute_velocity=True
                )

            res_r = rv_data[:, :3] - rv_oe[:, :3]

            if self.use_rsw:
                for ti in range(t_data.size):
                    M_i2rsw = self.rot_i2rsw(rv_oe[ti, :])
                    res_r[ti] = M_i2rsw @ res_r[ti]
        else:
            res_r = rv_data[:, :3]
            u_k = None

        # compute the Chebyshev coefficients
        if isinstance(L2_alpha, np.ndarray):
            ephem_out = np.zeros((L2_alpha.size, self.n_params))
            for i in range(L2_alpha.size):
                ephem_copy = deepcopy(ephem)
                try:
                    ephem_out[i] = self.poly_fit(
                        t_k,
                        res_r,
                        ephem_copy,
                        print_result,
                        u_k,
                        constrain_endpoints=False,
                        L2_alpha=L2_alpha[i],
                    )
                except Exception as e:
                    print(f"Error in fitting with L2_alpha={L2_alpha[i]}: {e}")
                    ephem_out[i] = np.nan * np.ones(self.n_params)
        else:
            ephem_out = self.poly_fit(
                t_k, res_r, ephem, print_result, u_k, constrain_endpoints=False, L2_alpha=L2_alpha
            )

        return ephem_out

    def plot_parmas(self, t_data, rvbf_data, ephem):
        """
        Plot the parameters of the ephemeris

        Args:
            t_data: array of times
            rvbf_data: array of position and velocity vectors in body frame
            rvbf_w_data: array of position and velocity vectors in body frame with velocity
        """

        t_ref, t_fit, ephem0 = self.init_guess(t_data, rvbf_data)

        # true params
        if self.use_meq:
            params_true = cart_to_mqoe(rvbf_data, self.GM)
        else:
            params_true = np.zeros((t_data.size, 4))
            for ti in range(t_data.size):
                coe = pnt.cart_to_classical(rvbf_data[ti], self.GM)
                a = coe[0]
                e = coe[1]
                inc = coe[2]
                lam = coe[3]
                w = coe[4]
                M = coe[5]
                nu = pnt.mean_to_true_anomaly(M, e)

                r = np.linalg.norm(rvbf_data[ti, :3])
                u = w + nu
                params_true[ti, :] = np.array([r, u, inc, lam])

        # compute params from the ephemeris
        res, params = self.ephem2cart(t_data, ephem, compute_velocity=True, return_params=True)

        # initialize the coefficients
        if self.use_meq:
            params_labels = ["p_k", "f_k", "g_k", "h_k", "k_k", "L_k"]
            rows = 2
            cols = 3
        else:
            params_labels = ["r_k", "u_k", "i_k", "lambda_k"]
            rows = 1
            cols = 4

        fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4))
        axes = axes.flatten()

        lw = 2

        for i in range(rows * cols):
            true_val = np.unwrap(params_true[:, i], np.pi)
            est_val = np.unwrap(params[params_labels[i]], np.pi)

            axes[i].plot(t_data, true_val, "r-", linewidth=lw, label="True")
            axes[i].plot(t_data, est_val, "b-", linewidth=lw, label="Fit")
            axes[i].legend()
            axes[i].set_xlabel("Time [min]")
            axes[i].set_ylabel(params_labels[i])
            axes[i].grid(True)

        plt.tight_layout()
        plt.show()
