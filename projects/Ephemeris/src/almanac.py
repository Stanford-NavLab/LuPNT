import numpy as np
import pylupnt as pnt
import matplotlib.pyplot as plt
from scipy.linalg import lstsq
import cvxpy as cp
from src.orbit_utils import convert_ci2pa, convert_pa2ci
from scipy.optimize import curve_fit

try:
    from .ephemeris import Ephemeris
    from .basis import get_poly_basis, get_fourier_basis, get_fourier_basis_dt, get_poly_basis_dt
    from .orbit_utils import wrapToPi, wrapTo2Pi, rms
except ImportError:
    from ephemeris import Ephemeris
    from basis import get_poly_basis, get_fourier_basis, get_fourier_basis_dt, get_poly_basis_dt
    from orbit_utils import wrapToPi, wrapTo2Pi, rms


def initialize_orbit(orbit, svid=0):

    if orbit == "LCRNS":
        # https://esc.gsfc.nasa.gov/static-files/LCRNS_Reference_Constellation_White_Paper_03_2025.pdf
        t0_tai = pnt.convert_time(pnt.gregorian_to_time(2027, 3, 1, 0, 0, 0), pnt.UTC, pnt.TAI)

        # a, e, i, Omega, w, M0
        svoe = np.zeros((5, 6))
        svoe[0] = np.array([11315.936501, 0.691982, 59.373229, 321.019197, 92.494031, 0.000000])
        svoe[1] = np.array([11317.948675, 0.691982, 58.951732, 320.997768, 92.505016, 180.000000])
        svoe[2] = np.array([11305.413654, 0.691982, 52.733096, 81.148790, 92.062891, 140.049207])
        svoe[3] = np.array([11326.302154, 0.691982, 52.513419, 81.138818, 92.068945, 195.992393])
        svoe[4] = np.array([11307.882863, 0.691982, 56.310396, 204.889626, 85.444071, 164.007607])
        svoe[:, 2:] = np.deg2rad(svoe[:, 2:])

        svoe[:, 0] *= 1e3  # km to m

        coe_pa = svoe[svid]
        a = coe_pa[0]
        ecc = coe_pa[1]
        inc = coe_pa[2]
        Omega = coe_pa[3]
        w = coe_pa[4]
        M0 = coe_pa[5]
        rv0_pa0 = pnt.classical_to_cart(coe_pa, pnt.GM_MOON)  # In PA frame
        rv0_ci = convert_pa2ci(t0_tai, rv0_pa0, rotate_only=True)

    elif orbit == "Moonlight":
        t0_tai = pnt.convert_time(pnt.gregorian_to_time(2027, 1, 1, 0, 0, 0), pnt.TDB, pnt.TAI)
        a = 9748.14e3  # meters
        ecc = 0.70
        inc = np.deg2rad(48.04)
        w = np.deg2rad(123.60)
        Omega = np.deg2rad(89.49)
        M0 = np.deg2rad(pnt.true_to_mean_anomaly(np.deg2rad(90.0), ecc))
        coe = np.array([a, ecc, inc, Omega, w, M0])

        rv0_ci = pnt.classical_to_cart(coe, pnt.GM_MOON)  # In CI frame

    elif orbit == "LNSS":
        t0_tai = pnt.convert_time(pnt.gregorian_to_time(2027, 1, 1, 0, 0, 0), pnt.TDB, pnt.TAI)
        a = 6541.4e3  # meters
        ecc = 0.60
        inc = np.deg2rad(56.2)
        w = np.deg2rad(90.0)
        Omega = np.deg2rad(0.0)
        M0 = np.deg2rad(pnt.true_to_mean_anomaly(np.deg2rad(0.0), ecc))
        coe_op = np.array([a, ecc, inc, Omega, w, M0])
        rv0_op = pnt.classical_to_cart(coe_op, pnt.GM_MOON)
        rv0_ci = pnt.convert_frame(t0_tai, rv0_op, pnt.MOON_OP, pnt.MOON_CI)  # In CI frame

    elif orbit == "Polar":
        t0_tai = pnt.convert_time(pnt.gregorian_to_time(2027, 1, 1, 0, 0, 0), pnt.TDB, pnt.TAI)

        a = 3870.0e3
        ecc = 0.0
        inc = np.deg2rad(104.428)
        Omega = np.deg2rad(53.563)
        w = np.deg2rad(90)
        M0 = np.deg2rad(pnt.true_to_mean_anomaly(-5.0, ecc))

        coe_op = np.array([a, ecc, inc, Omega, w, M0])
        rv0_ci = pnt.classical_to_cart(coe_op, pnt.GM_MOON)
        # rv0_ci = pnt.convert_frame(t0_tai, rv0_op, pnt.MOON_OP, pnt.MOON_CI)  # In CI frame

    elif orbit == "NRHO":
        t0_tai = pnt.convert_time(pnt.gregorian_to_time(2022, 7, 15, 0, 0, 0), pnt.UTC, pnt.TAI)
        rv0_ci = (
            np.array(
                [
                    4811.24915,
                    23140.1109,
                    -66343.1582,
                    -0.0560634369,
                    -0.0401916424,
                    -0.0180843762,
                ]
            )
            * 1e3
        )
        T_syn = 29.68278536728694 * 86400
        # synodic period
        period = T_syn * 2 / 9
        twobody = False

        coe = pnt.cart_to_classical(rv0_ci, pnt.GM_MOON)
        a = coe[0]
        ecc = coe[1]
        inc = coe[2]
    print("Initial position and velocity in CI frame:")
    print(rv0_ci)

    return t0_tai, rv0_ci


def setup_orbit(t0_tai, rv0_mci, t_days, plot_cart=False, plot_coe=False, figname=None):

    # dynamics
    add_earth = True
    add_sun = True
    sphm = [20, 20]  # Spherical harmonic model degree and order

    # dynamics
    dyn = pnt.NBodyDynamics()
    dyn.set_integrator(pnt.IntegratorType.RKF45)
    dyn.set_integrator_params(pnt.IntegratorParams(max_iter=20, abstol=1e-12, reltol=1e-12))
    dyn.add_body(pnt.Body.Moon(sphm[0], sphm[1]))
    if add_earth:
        dyn.add_body(pnt.Body.Earth())
    if add_sun:
        dyn.add_body(pnt.Body.Sun())
    dyn.set_frame(pnt.MOON_CI)
    dyn.set_time_step(30.0)  # propagation timestep

    # propagate
    print("Propagating orbit for {} days ...".format(t_days))
    tend = t_days * 86400.0  # days to seconds
    tspan = np.linspace(0.0, tend, int(tend / 600) + 1)  # every minute
    t_tai = tspan + t0_tai
    pnt.set_lupnt_epoch(0)
    rv_prop_mci = dyn.propagate(rv0_mci, t_tai)
    coe_prop_mci = pnt.cart_to_classical(rv_prop_mci, pnt.GM_MOON)

    # convert to PA frame
    print("Converting to PA frame ...")
    t_tai = tspan + t0_tai
    rv_prop_pa = convert_ci2pa(t_tai, rv_prop_mci, rotate_only=True)
    rv_prop_pa_w = convert_ci2pa(t_tai, rv_prop_mci, rotate_only=False)
    coe_prop_pa = pnt.cart_to_classical(rv_prop_pa, pnt.GM_MOON)  # In PA frame

    rv_prop_op = pnt.convert_frame(t_tai, rv_prop_mci, pnt.MOON_CI, pnt.MOON_OP)
    coe_prop_op = pnt.cart_to_classical(rv_prop_op, pnt.GM_MOON)  # In OP frame

    # convert to classical orbital elements
    if plot_cart:
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
        labels = ["x", "y", "z", "vx", "vy", "vz"]
        for i in range(3):
            axes[0, i].plot((t_tai - t0_tai) / 86400, rv_prop_mci[:, i] / 1e3)
            axes[0, i].set_ylabel(f"{labels[i]} [km]")
            axes[0, i].grid(True)

            axes[1, i].plot((t_tai - t0_tai) / 86400, rv_prop_mci[:, i + 3] / 1e3)
            axes[1, i].set_ylabel(f"{labels[i+3]} [km/s]")
            axes[1, i].grid(True)

            axes[1, i].set_xlabel("Days")

        plt.tight_layout()

    if plot_coe:
        fig, axes = plt.subplots(1, 5, figsize=(20, 4))
        # plot_list = [coe_prop_mci, coe_prop_pa]
        # labels = ["CI", "PA"]
        plot_list = [coe_prop_pa, coe_prop_op]
        labels = ["PA", "OP"]

        for j, coe_prop in enumerate(plot_list):
            label = labels[j]
            coe_plot = np.zeros_like(coe_prop)
            coe_plot[:, 0] = coe_prop[:, 0] / 1e3  # a in km, e unitless
            coe_plot[:, 1] = coe_prop[:, 1]  # e
            coe_plot[:, 2:] = np.rad2deg(coe_prop[:, 2:])

            labels = [
                "Semi-major Axis [km]",
                "Eccentricity [] ",
                "Inclination [deg]",
                "Longitude of Ascending Node [deg]",
                "Argument of Periapsis [deg]",
                "Mean Anomaly [deg]",
            ]
            for i in range(5):
                axes[i].plot((t_tai - t0_tai) / 86400, coe_plot[:, i], label=label)
                axes[i].set_ylabel(f"{labels[i]}", fontsize=14)
                axes[i].grid(True)
                axes[i].set_xlabel("Days", fontsize=14)
            # axes[0, i].legend()
            # axes[1, i].legend()
        plt.tight_layout()

        if figname is not None:
            fig.savefig(figname, dpi=300)

    results = {
        "t_tai": t_tai,
        "rv_prop_mci": rv_prop_mci,
        "rv_prop_pa": rv_prop_pa,
        "rv_prop_pa_w": rv_prop_pa_w,
        "coe_prop_mci": coe_prop_mci,
        "coe_prop_pa": coe_prop_pa,
        "coe_prop_op": coe_prop_op,
    }

    return results


class Almanac(Ephemeris):
    def __init__(
        self,
        param_dict,
        body=pnt.MOON,
        print_info=False,
        use_sma_special=False,
        use_full_sidereal=False,
        use_op=False,
    ):
        super().__init__(body)
        self.param_dict = param_dict
        self.use_sma_special = use_sma_special

        self.create_ephem_dict(print_info=print_info)
        # sidereal period
        self.omega_sidereal = 2 * np.pi / self.T_sidereal
        if use_full_sidereal:
            self.omega_sidereal = np.pi / self.T_sidereal_full
        self.use_op = use_op

    def create_ephem_dict(self, print_info=False):

        self.params = ["a", "e", "w", "M", "i", "l", "u"]
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

                if self.use_sma_special and param_str == "a":
                    print("Using special SMA fitting ...")
                    keys_sma = ["a_0", "a_1"]
                    for i in range(self.poly_orders["a"]):
                        for j in range(self.fourier_orders["a"]):
                            keys_f = [f"a_{i}_C{j}", f"a_{i}_S{j}"]
                            keys_sma.extend(keys_f)
                    for key in keys_sma:
                        self.idx_keys[key] = len(self.keys_list)
                        self.keys_list.append(key)

                else:
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

    def get_sma_design_matrix(self, t_k, Phi, dict):
        param_str = "a"
        n = t_k.size
        Amat_ps = get_fourier_basis(Phi, self.fourier_orders[param_str])  # N x (2*fsn)
        Amat_poly = get_poly_basis(
            t_k, self.poly_orders[param_str], dict["t_fit"], self.poly_types[param_str]
        )  # N x (ln+1)
        # y = a0 + a1*t + (b0 + b1*t + ...) * (C1 cos(w_sid t) + S1 sin(w_sid t) + ...)
        # Expand to get the design matrix
        n_poly = Amat_poly.shape[1]
        n_fourier = int(Amat_ps.shape[1] / 2)
        Amat = np.zeros((n, n_poly * n_fourier * 2 + 2))  # +2 for a0 and a1
        keys = []
        # constant term a0
        Amat[:, 0] = 1.0
        keys.append("a_0")
        # linear term a1
        Amat[:, 1] = t_k / dict["t_fit"]
        keys.append("a_1")
        # cross terms
        col_idx = 2
        for i in range(n_poly):
            for j in range(n_fourier):
                Amat[:, col_idx] = Amat_poly[:, i] * Amat_ps[:, j]
                keys_f = [f"a_{i}_C{j}", f"a_{i}_S{j}"]
                keys.extend(keys_f)
                col_idx += 1

        return Amat, keys

    def compute_param(self, dict, param_str, t_k, Phi, offset):
        """
        Compute the value of the parameter at time t_k
        """
        pvalue = offset

        if self.use_sma_special and param_str == "a":
            T, keys = self.get_sma_design_matrix(t_k, Phi, dict)
            a_coeffs = np.array([self.get_value(dict, key, 0) for key in keys])
            pvalue += T @ a_coeffs
        else:
            if self.use_poly:
                Tp = get_poly_basis(t_k, self.max_poly_order, dict["t_fit"], "monomial")
            if self.use_cheby:
                Tc = get_poly_basis(t_k, self.max_poly_order, dict["t_fit"], "chebyshev")
            if self.use_fourier and Phi is not None:
                Tf = get_fourier_basis(Phi, self.max_fourier_order)
            if self.use_fourier_sidereal:
                Tfs = get_fourier_basis(self.omega_sidereal * t_k, self.max_fourier_sidereal_order)

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

        # pre-computation for chebyshev
        if self.use_sma_special and param_str == "a":
            T, keys = self.get_sma_design_matrix(t_k, Phi, dict)
            a_coeffs = np.array([self.get_value(dict, key, 0) for key in keys])
            pvalue += T @ a_coeffs
        else:
            if self.use_poly:
                Tpdot = get_poly_basis_dt(t_k, self.max_poly_order, dict["t_fit"], "monomial")
            if self.use_cheby:
                Tcdot = get_poly_basis_dt(t_k, self.max_poly_order, dict["t_fit"], "chebyshev")
            if self.use_fourier and Phi is not None:
                Phidot = self.compute_Phidot(dict, t_k)
                Tfdot = get_fourier_basis_dt(Phi, Phidot, self.max_fourier_order)
            if self.use_fourier_sidereal:
                Tfsdot = get_fourier_basis_dt(
                    self.omega_sidereal * t_k, self.omega_sidereal, self.max_fourier_sidereal_order
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
        """
        Compute the time derivative of Phi
        """
        a_0 = dict["a_0"]
        T_orbit = 4 * np.pi * np.sqrt(a_0**3 / self.GM)
        Phidot = 2 * np.pi / T_orbit * np.ones_like(t_k)
        return Phidot

    def compute_rdot(self, dict, t_k):
        """
        Compute the radius at time t_k
        """
        a_0 = dict["a_0"]
        T_orbit = 4 * np.pi * np.sqrt(a_0**3 / self.GM)
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

    def init_guess(self, t_data, rvbf_data, print_result=False):
        """
        Initialize the coefficients for the alamanc

        Args:
            t_data: array of times
            rv_data: array of position and velocity vectors in body frame
        """
        lent = t_data.size
        idx = 0  # use the initial index
        t_ref = t_data[idx]
        rv_ref = rvbf_data[idx]

        t_fit = t_data[-1] - t_data[0]  # time span

        if print_result:
            print("  [Initial guess] ")
            print("    t_ref: ", t_ref)
            print("    t_fit: ", t_fit)

        if self.use_op:
            rv_ref = pnt.convert_frame(t_ref, rv_ref, pnt.MOON_PA, pnt.MOON_OP, rotate_only=True)

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
            "i": np.zeros(t_data.size),
            "l": np.zeros(t_data.size),
            "w": np.zeros(t_data.size),
            "M": np.zeros(t_data.size),
            "u": np.zeros(t_data.size),
            "f": np.zeros(t_data.size),
        }
        print("shapes: ", rv_data.shape, t_data.shape)

        if self.use_op:
            for i, t in enumerate(t_data):
                rv_data[i] = pnt.convert_frame(
                    t, rv_data[i], pnt.MOON_PA, pnt.MOON_OP, rotate_only=True
                )

        coe = pnt.cart_to_classical(rv_data, self.GM)  # N x 6
        a = coe[:, 0]
        e = coe[:, 1]
        inc = coe[:, 2]
        lamda = wrapTo2Pi(coe[:, 3])
        w = coe[:, 4]
        M = wrapTo2Pi(coe[:, 5])
        f = wrapTo2Pi(pnt.mean_to_true_anomaly(coe[:, 5], coe[:, 1]))
        u = wrapTo2Pi(w + f)

        # make sure the angles are continuous
        M = np.unwrap(M)
        w = np.unwrap(w)
        inc = np.unwrap(inc)
        lamda = np.unwrap(lamda)

        param_data["a"] = a
        param_data["e"] = e
        param_data["i"] = inc
        param_data["l"] = lamda
        param_data["w"] = w
        param_data["M"] = M
        param_data["u"] = u
        param_data["f"] = f

        print("shapes after: ", param_data["a"].shape, param_data["M"].shape)

        return param_data

    def fit_coeffs(self, residuals, param_str, ephem, t_k, Phi, print_result=False, plot_fit=False):
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

        # setup the design matrix ---------------------------------
        if self.use_sma_special and param_str == "a":
            Amat, keys = self.get_sma_design_matrix(t_k, Phi, dict)

        else:
            # polynomials
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
                Amat_ps = get_fourier_basis(
                    self.omega_sidereal * t_k, self.fourier_sidereal_orders[param_str]
                )
                Amat = np.hstack((Amat, Amat_ps))
                keys.extend(self.get_fourier_sidereal_coeff_keys(param_str))

        # Solve Optimization problem ---------------------------------
        n = Amat.shape[0]

        # solve using scipy least squares
        if self.fit_obj == "lsq":
            coeffs, res_opt, rank, s = lstsq(Amat, residuals)
        else:
            n = Amat.shape[1]
            x = cp.Variable(n)

            if self.fit_obj == "max":
                obj = cp.Minimize(cp.norm(Amat @ x - residuals, "inf"))
            elif self.fit_obj == "lsq-cvx":
                obj = cp.Minimize(cp.norm(Amat @ x - residuals, 2))
            else:
                raise ValueError("Invalid fit objective: {}".format(self.fit_obj))

            prob = cp.Problem(obj, [])
            prob.solve(solver=cp.ECOS, max_iters=10000)
            coeffs = x.value
            res_opt = np.sum((Amat @ coeffs - residuals) ** 2)

        res_rms = np.sqrt(res_opt / n)
        if param_str == "M":
            res_rms = np.mod(res_rms, np.pi * 2)
        if print_result:
            print("  [Fit {}] ".format(param_str))
            print("    Residual after fitting: ", res_rms)
            print("    Coefficients: ", coeffs)

        # replace the coefficients in the ephemeris -----------------------
        ephem = self.replace_ephem(ephem, keys, coeffs)

        # plot the residuals and fit
        plot_diff = False  # plot data - fit

        param_str_long = {
            "a": "Semi-major Axis (km)",
            "e": "Eccentricity",
            "i": "Inclination (rad)",
            "l": "Longitude of Ascending Node (rad)",
            "w": "Argument of Periapsis (rad)",
            "M": "Mean Anomaly (rad)",
            "u": "Argument of Latitude (rad)",
        }.get(param_str, param_str)

        if plot_fit:
            fitvals = Amat @ coeffs
            plt.figure(figsize=(3, 3))
            if param_str == "a":
                residuals = (residuals - coeffs[0]) / 1000
                fitvals = (fitvals - coeffs[0]) / 1000
            if param_str == "M":
                residuals = wrapToPi(residuals)
                fitvals = wrapToPi(fitvals)

            if plot_diff:
                plt.plot(
                    (t_k + dict["t_ref"]) / 86400,
                    residuals - fitvals,
                    "k.",
                    label="data - fit",
                    markersize=4,
                )
            else:
                plt.plot((t_k + dict["t_ref"]) / 86400, residuals, "k.", label="data", markersize=4)
                plt.plot((t_k + dict["t_ref"]) / 86400, fitvals, "r-", label="fit", linewidth=2)
            plt.xlabel("Time (days)")
            plt.ylabel("{} residual".format(param_str))
            plt.title("{}".format(param_str_long))

            # for k in range(1):
            #     plt.axvline(x=(k+1) * self.T_sidereal/2 / 86400, color="gray", linestyle="--")

            # T_orbit = 2 * np.pi * np.sqrt(self.get_value(self.ephem2dict(ephem), "a_0", 1.0e7)**3 / self.GM)
            # for k in range(1):
            #     plt.axvline(x=(k+1) * 2 * T_orbit / 86400,  color="blue", linestyle="--")

            plt.legend()
            plt.grid()
            plt.tight_layout()
            if self.figname is not None:
                plt.savefig(f"{self.figname}_{param_str}_fit.pdf")
            plt.show()

        return ephem

    def _cumtrapz(self, y, t):
        """Simple cumulative trapezoid integral with y(t) samples."""
        y = np.asarray(y, float)
        t = np.asarray(t, float)
        out = np.zeros_like(y)
        dt = np.diff(t)
        # trapezoid for each step k: 0.5*(y[k]+y[k-1])*(t[k]-t[k-1])
        incr = 0.5 * (y[1:] + y[:-1]) * dt
        out[1:] = np.cumsum(incr)
        return out  # integral from t[0] to each t[k]

    def fit_almanac_coeffs(self, t_data, param_data, almanac, print_result=False, plot_fit=False):
        dict = self.ephem2dict(almanac)

        # extract the fixed parameters
        t_ref = dict["t_ref"]
        t_k = t_data - t_ref
        ndata = t_k.size

        a_0 = dict["a_0"]
        T_orbit = 4 * np.pi * np.sqrt(a_0**3 / self.GM)
        Phi = 2 * np.pi / T_orbit * t_k  # tempoaral variable for fourier

        # 1. fit semi-major axis
        almanac = self.fit_coeffs(
            param_data["a"], "a", almanac, t_k, Phi, print_result=print_result, plot_fit=plot_fit
        )
        a_k = self.compute_param(self.ephem2dict(almanac), "a", t_k, Phi, 0)

        # update a_0
        a_0 = self.get_value(self.ephem2dict(almanac), "a_0", a_0)
        T_orbit = 4 * np.pi * np.sqrt(a_0**3 / self.GM)
        Phi = 2 * np.pi / T_orbit * t_k  # tempoaral variable for fourier

        # 2. fit eccentricity
        almanac = self.fit_coeffs(
            param_data["e"], "e", almanac, t_k, Phi, print_result=print_result, plot_fit=plot_fit
        )
        e_k = self.compute_param(self.ephem2dict(almanac), "e", t_k, Phi, 0)

        # 3. fit inclination
        alamanc = self.fit_coeffs(
            param_data["i"], "i", almanac, t_k, Phi, print_result=print_result, plot_fit=plot_fit
        )

        # 4. Argument of periapsis
        # almanac = self.fit_coeffs(
        #     param_data["w"], "w", almanac, t_k, Phi, print_result=print_result, plot_fit=plot_fit
        # )
        # w_k = self.compute_param(self.ephem2dict(almanac), "w", t_k, Phi, 0)
        w_k = 0.0

        # 5. Longitude of ascending node
        res_lambda = np.unwrap(param_data["l"] + self.omega_b * t_k)
        almanac = self.fit_coeffs(
            res_lambda, "l", almanac, t_k, Phi, print_result=print_result, plot_fit=plot_fit
        )

        # 6. Mean anomaly
        n_t = np.sqrt(self.GM * np.ones_like(a_k) / a_k**3)
        M_nom = self._cumtrapz(n_t, t_k)

        res_M = np.unwrap(param_data["M"] - M_nom)
        almanac = self.fit_coeffs(
            res_M, "M", almanac, t_k, Phi, print_result=print_result, plot_fit=plot_fit
        )
        M_k = self.compute_param(self.ephem2dict(almanac), "M", t_k, Phi, M_nom)

        E_k = M_k
        for _ in range(5):
            E_k = E_k - (E_k - e_k * np.sin(E_k) - M_k) / (1 - e_k * np.cos(E_k))
        nu_k = 2 * np.arctan(np.sqrt((1 + e_k) / (1 - e_k)) * np.tan(E_k / 2))  # true anomaly

        # 7. argument of latitude
        res_u = np.unwrap(param_data["u"] - w_k - nu_k)
        almanac = self.fit_coeffs(
            res_u, "u", almanac, t_k, Phi, print_result=print_result, plot_fit=plot_fit
        )

        return almanac

    def fit(self, t_data, rv_data, fit_obj="lsq", print_result=False, plot_fit=False, figname=None):

        self.fit_obj = fit_obj
        self.figname = figname

        t_ref, t_fit, almanac = self.init_guess(t_data, rv_data, print_result=print_result)

        param_data = self.generate_param_data(t_data, rv_data)

        alamanc = self.fit_almanac_coeffs(
            t_data, param_data, almanac, print_result=print_result, plot_fit=plot_fit
        )

        return alamanc

    def ephem2cart(self, t, ephem, compute_velocity=False, return_params=False):

        # base params --------------------------------
        dict = self.ephem2dict(ephem)
        t_ref = dict["t_ref"]
        t_k = t - t_ref
        a_0 = dict["a_0"]
        T_orbit = 4 * np.pi * np.sqrt(a_0**3 / self.GM)
        Phi = 2 * np.pi / T_orbit * t_k  # tempoaral variable for fourier
        n = np.sqrt(self.GM / a_0**3)

        # orbital elements --------------------------------------------------
        a_k = self.compute_param(dict, "a", t_k, Phi, 0)
        e_k = self.compute_param(dict, "e", t_k, Phi, 0)  # eccentricity
        i_k = self.compute_param(dict, "i", t_k, Phi, 0)  # inclination
        lambda_k = self.compute_param(
            dict, "l", t_k, Phi, -self.omega_b * t_k
        )  # longitude of ascending node
        # w_k = self.compute_param(dict, "w", t_k, Phi, 0)  # argument of periapsis
        w_k = 0.0

        n_t = np.sqrt(self.GM * np.ones_like(a_k) / a_k**3)
        M_nom = self._cumtrapz(n_t, t_k)
        M_k = self.compute_param(dict, "M", t_k, Phi, M_nom)  # mean anomaly

        # solve Kepler's equation E = M + e * sin(E) -----------------------------
        E_k = M_k
        for _ in range(5):
            E_k = E_k - (E_k - e_k * np.sin(E_k) - M_k) / (1 - e_k * np.cos(E_k))
        nu_k = 2 * np.arctan(np.sqrt((1 + e_k) / (1 - e_k)) * np.tan(E_k / 2))  # true anomaly

        r_k = a_k * (1 - e_k * np.cos(E_k))  # radius
        u_k = w_k + nu_k

        u_k = self.compute_param(dict, "u", t_k, Phi, u_k)  # argument of latitude

        # positions --------------------------------------------------
        x_k = r_k * np.cos(u_k)
        y_k = r_k * np.sin(u_k)
        x_pos = x_k * np.cos(lambda_k) - y_k * np.cos(i_k) * np.sin(lambda_k)
        y_pos = x_k * np.sin(lambda_k) + y_k * np.cos(i_k) * np.cos(lambda_k)
        z_pos = y_k * np.sin(i_k)

        # velocities --------------------------------------------------
        if compute_velocity:

            # for r and u, use diff
            rdot0 = self.compute_rdot(dict, t_k)
            Phidot = self.compute_Phidot(dict, t_k)
            rdot = rdot0
            udot = Phidot
            didt = self.compute_param_dt(dict, "i", t_k, Phi, 0)  # inclination
            lambdadot = self.compute_param_dt(
                dict, "l", t_k, Phi, -self.omega_b
            )  # longitude of ascending node

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

        if self.use_op:
            res = pnt.convert_frame(t, res, pnt.MOON_OP, pnt.MOON_PA, rotate_only=False)

        return res
