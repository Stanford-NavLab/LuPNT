import numpy as np
from scipy.linalg import lstsq

try:
    from .ephemeris import Ephemeris
except ImportError:
    from ephemeris import Ephemeris

import pylupnt as pnt
from scipy.optimize import least_squares


class GPSEphemeris(Ephemeris):
    """
    GPS/Galileo broadcast ephemeris coordinate computation (ECEF/TRS) per Navipedia.

    Parameters included (position):
      toe, sqrtA, e, M0, omega, i0, Omega0, DeltaN, iDot, OmegaDot,
      Cuc, Cus, Crc, Crs, Cic, Cis

    Notes:
      - Use body=pnt.EARTH so GM and omega_b match Earth.
      - Returns ECEF at the *transmission-time* frame (as in the Navipedia algorithm).
        If you need "ECEF at reception time", you must apply the additional Earth-rotation
        correction separately (Navipedia has a dedicated article for that).
    """

    def __init__(self, body, print_info=False):
        super().__init__(body=body)
        self.create_ephem_dict(print_info=print_info)

    def create_ephem_dict(self, print_info=False):
        self.idx_keys = {}
        self.keys_list = []

        # Keep a t_ref slot for compatibility with your other ephemeris families.
        # Here we don't use it inside ephem2cart; toe is the real broadcast reference.
        self.idx_keys["t_ref"] = len(self.keys_list)
        self.keys_list.append("t_ref")

        params = [
            "sqrtA",
            "e",
            "M0",
            "omega",
            "i0",
            "Omega0",
            "DeltaN",
            "iDot",
            "OmegaDot",
            "Cuc",
            "Cus",
            "Crc",
            "Crs",
            "Cic",
            "Cis",
        ]
        for p in params:
            self.idx_keys[p] = len(self.keys_list)
            self.keys_list.append(p)

        self.n_params = len(self.keys_list)

        if print_info:
            print("key list:", self.keys_list)
            print("keys idx :", self.idx_keys)
            print("Number of parameters:", self.n_params)

    @staticmethod
    def _week_crossover(tk):
        """
        GPS week crossover adjustment:
          if tk > 302400 -> tk -= 604800
          if tk < -302400 -> tk += 604800
        """
        tk = np.asarray(tk, dtype=float)
        tk = np.where(tk > 302400.0, tk - 604800.0, tk)
        tk = np.where(tk < -302400.0, tk + 604800.0, tk)
        return tk

    @staticmethod
    def _solve_kepler(M, e, max_iter=20, tol=1e-13):
        """
        Solve M = E - e*sin(E) for E using Newton iterations.
        Works with scalar or ndarray.
        """
        M = np.asarray(M, dtype=float)
        E = M.copy()
        for _ in range(max_iter):
            f = E - e * np.sin(E) - M
            fp = 1.0 - e * np.cos(E)
            dE = -f / fp
            E = E + dE
            if np.max(np.abs(dE)) < tol:
                break
        return E

    def ephem2cart(self, t, ephem, compute_velocity=False, return_params=False, scale=None):
        """
        Compute satellite ECEF position from broadcast ephemeris at epochs t.

        Args:
          t: scalar or (N,) time [s] in GPS/Galileo system time consistent with toe
          ephem: (n_params,) array
          compute_velocity: if True, appends finite-difference velocity [m/s]
          return_params: if True, also return intermediate variables dict
          scale: optional quantization scale (passed through ephem2dict)

        Returns:
          pos: (N,3) or (N,6) array
          (optional) params dict
        """
        d = self.ephem2dict(ephem, scale=scale)

        t = np.atleast_1d(np.asarray(t, dtype=float))
        toe = float(d["t_ref"])
        tk = self._week_crossover(t - toe)

        # Broadcast params
        sqrtA = float(d["sqrtA"])
        e = float(d["e"])
        M0 = float(d["M0"])
        omega = float(d["omega"])
        i0 = float(d["i0"])
        Omega0 = float(d["Omega0"])

        DeltaN = float(d["DeltaN"])
        iDot = float(d["iDot"])
        OmegaDot = float(d["OmegaDot"])

        Cuc = float(d["Cuc"])
        Cus = float(d["Cus"])
        Crc = float(d["Crc"])
        Crs = float(d["Crs"])
        Cic = float(d["Cic"])
        Cis = float(d["Cis"])

        # Constants
        mu = float(self.GM)
        omega_E = float(self.omega_b)

        # a, mean motion
        a = sqrtA * sqrtA
        n0 = np.sqrt(mu / (a * a * a))
        n = n0 + DeltaN

        # Mean anomaly
        Mk = M0 + n * tk

        # Eccentric anomaly
        Ek = self._solve_kepler(Mk, e)

        # True anomaly (atan2 form)
        sin_vk = np.sqrt(1.0 - e * e) * np.sin(Ek) / (1.0 - e * np.cos(Ek))
        cos_vk = (np.cos(Ek) - e) / (1.0 - e * np.cos(Ek))
        vk = np.arctan2(sin_vk, cos_vk)

        # Argument of latitude, corrections
        phi = omega + vk
        two_phi = 2.0 * phi

        du = Cuc * np.cos(two_phi) + Cus * np.sin(two_phi)
        dr = Crc * np.cos(two_phi) + Crs * np.sin(two_phi)
        di = Cic * np.cos(two_phi) + Cis * np.sin(two_phi)

        uk = phi + du
        rk = a * (1.0 - e * np.cos(Ek)) + dr
        ik = i0 + iDot * tk + di

        # Longitude of ascending node (Greenwich)
        lambdak = Omega0 + (OmegaDot - omega_E) * tk - omega_E * toe

        # Orbital plane coords
        x_orb = rk * np.cos(uk)
        y_orb = rk * np.sin(uk)

        cosO = np.cos(lambdak)
        sinO = np.sin(lambdak)
        cosi = np.cos(ik)
        sini = np.sin(ik)

        # ECEF
        x = x_orb * cosO - y_orb * cosi * sinO
        y = x_orb * sinO + y_orb * cosi * cosO
        z = y_orb * sini

        pos = np.vstack([x, y, z]).T  # (N,3)

        if not compute_velocity:
            if return_params:
                params = dict(
                    tk=tk, Mk=Mk, Ek=Ek, vk=vk, phi=phi, uk=uk, rk=rk, ik=ik, lambdak=lambdak
                )
                return pos, params
            return pos

        # Finite-difference velocity (simple + robust)
        dt = 0.1  # seconds; tune if needed
        pos_p = self.ephem2cart(
            t + dt, ephem, compute_velocity=False, return_params=False, scale=scale
        )
        pos_m = self.ephem2cart(
            t - dt, ephem, compute_velocity=False, return_params=False, scale=scale
        )
        vel = (pos_p - pos_m) / (2.0 * dt)

        out = np.hstack([pos, vel])  # (N,6)
        if return_params:
            params = dict(tk=tk, Mk=Mk, Ek=Ek, vk=vk, phi=phi, uk=uk, rk=rk, ik=ik, lambdak=lambdak)
            return out, params
        return out

    def coe2dict(self, coe, t_ref):
        dict = {}

        keys_to_fit = [
            "t_ref",
            "sqrtA",
            "e",
            "M0",
            "omega",
            "i0",
            "Omega0",
            "DeltaN",
            "iDot",
            "OmegaDot",
            "Cuc",
            "Cus",
            "Crc",
            "Crs",
            "Cic",
            "Cis",
        ]

        for i, key in enumerate(self.keys_list):
            if key == "t_ref":
                dict[key] = t_ref
            elif key == "sqrtA":
                dict[key] = np.sqrt(coe[0])
            elif key == "e":
                dict[key] = coe[1]
            elif key == "i0":
                dict[key] = coe[2]
            elif key == "Omega0":
                dict[key] = coe[3]
            elif key == "omega":
                dict[key] = coe[4]
            elif key == "M0":
                dict[key] = coe[5]
            else:
                dict[key] = 0

        return dict

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

        if print_result:
            print("  [Initial guess] ")
            print("    t_ref: ", t_ref)

        coe = pnt.cart_to_classical(rv_ref, self.GM)
        dict = self.coe2dict(coe, t_ref)

        if print_result:
            print(dict)

        # initialize the coefficients
        ephem = np.array(self.dict2ephem(dict))

        return t_ref, ephem

    def fit_coeff(
        self,
        t_data,
        rv_ecef_data,
        max_iter=15,
        tol=1e-6,
        fd_rel_step=1e-7,
        fd_abs_floor=1e-10,
        damping=0.0,
        verbose=False,
    ):
        """
        Gauss–Newton fit of broadcast parameters to ECEF position samples.

        Args:
          t_data: (N,) times [s]
          r_ecef_data: (N,3) ECEF positions [m]
          ephem0: (n_params,) initial ephemeris vector
          keys_to_fit: list of ephemeris keys to estimate; default fits all broadcast keys (except t_ref)
          damping: Levenberg diagonal damping >=0 (0 disables)
        Returns:
          ephem: fitted ephemeris vector (same length as ephem0)
        """
        t_data = np.asarray(t_data, dtype=float).ravel()
        r_ecef_data = rv_ecef_data[:, :3]  # Use position part for fitting
        if r_ecef_data.shape != (t_data.size, 3):
            raise ValueError("r_ecef_data must have shape (N,3) matching t_data")

        keys_to_fit = [
            "t_ref",
            "sqrtA",
            "e",
            "M0",
            "omega",
            "i0",
            "Omega0",
            "DeltaN",
            "iDot",
            "OmegaDot",
            "Cuc",
            "Cus",
            "Crc",
            "Crs",
            "Cic",
            "Cis",
        ]

        fit_indices = [self.get_index(k) for k in keys_to_fit]
        if any(idx is None for idx in fit_indices):
            raise ValueError("One or more keys_to_fit not found in idx_keys")

        t_ref, ephem = self.init_guess(t_data, rv_ecef_data, verbose)

        def residual_vec(eph):
            pred = self.ephem2cart(t_data, eph, compute_velocity=False, return_params=False)
            return (pred - r_ecef_data).reshape(-1)  # (3N,)

        for it in range(max_iter):
            r = residual_vec(ephem)
            J = np.zeros((r.size, len(fit_indices)), dtype=float)

            # Finite-difference Jacobian
            for j, idx in enumerate(fit_indices):
                p0 = ephem[idx]
                step = max(fd_abs_floor, fd_rel_step * max(1.0, abs(p0)))
                ep_p = ephem.copy()
                ep_m = ephem.copy()
                ep_p[idx] = p0 + step
                ep_m[idx] = p0 - step

                rp = residual_vec(ep_p)
                rm = residual_vec(ep_m)
                J[:, j] = (rp - rm) / (2.0 * step)

            # Solve GN step
            if damping > 0.0:
                JTJ = J.T @ J
                g = J.T @ r
                dp = np.linalg.solve(JTJ + damping * np.eye(JTJ.shape[0]), -g)
            else:
                dp, *_ = lstsq(J, -r)

            # Update
            for j, idx in enumerate(fit_indices):
                ephem[idx] += dp[j]

            if verbose:
                rms = np.sqrt(np.mean(r * r))
                print(f"[iter {it:02d}] rms={rms:.6e}  ||dp||_inf={np.max(np.abs(dp)):.3e}")

            if np.max(np.abs(dp)) < tol:
                break

        return ephem

    def fit_coeff_scipy(
        self,
        t_data,
        rv_ecef_data,
        loss="soft_l1",
        f_scale=100.0,
        method="trf",
        max_nfev=200,
        x_scale="jac",
        verbose=2,
        bounds_dict=None,
        robust=True,
    ):
        """
        Robust nonlinear least squares fit using scipy.optimize.least_squares.

        Args:
          t_data: (N,) [s]
          rv_ecef_data: (N,6) [m, m/s]
          ephem0: (n_params,) initial guess
          keys_to_fit: list of keys to estimate; default fits all broadcast keys except 't_ref'
          loss: 'linear','soft_l1','huber','cauchy','arctan' (robust if not 'linear')  [oai_citation:1‡docs.scipy.org](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html?utm_source=chatgpt.com)
          f_scale: scale separating inliers/outliers for robust loss (units of residuals; meters here)  [oai_citation:2‡docs.scipy.org](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html?utm_source=chatgpt.com)
          method: 'trf' (recommended), 'dogbox', or 'lm' (lm only supports linear loss)  [oai_citation:3‡docs.scipy.org](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html?utm_source=chatgpt.com)
          bounds_dict: optional dict {key: (lb, ub)} in physical units; any missing keys are unbounded
          x_scale: 'jac' usually good when parameters have mixed units
        """
        t_data = np.asarray(t_data, dtype=float).ravel()
        rv_ecef_data = np.asarray(rv_ecef_data, dtype=float)
        if rv_ecef_data.shape != (t_data.size, 6):
            raise ValueError("rv_ecef_data must have shape (N,6) matching t_data")

        keys_to_fit = [
            "t_ref",
            "sqrtA",
            "e",
            "M0",
            "omega",
            "i0",
            "Omega0",
            "DeltaN",
            "iDot",
            "OmegaDot",
            "Cuc",
            "Cus",
            "Crc",
            "Crs",
            "Cic",
            "Cis",
        ]

        fit_indices = [self.get_index(k) for k in keys_to_fit]
        if any(idx is None for idx in fit_indices):
            raise ValueError("One or more keys_to_fit not found in idx_keys")

        t_ref, ephem0 = self.init_guess(t_data, rv_ecef_data, verbose)
        x0 = ephem0[fit_indices].copy()

        # bounds
        if bounds_dict is None:
            lb = -np.inf * np.ones_like(x0)
            ub = np.inf * np.ones_like(x0)
        else:
            lb = -np.inf * np.ones_like(x0)
            ub = np.inf * np.ones_like(x0)
            for j, k in enumerate(keys_to_fit):
                if k in bounds_dict and bounds_dict[k] is not None:
                    lb[j], ub[j] = bounds_dict[k]

        # If user requests robust but accidentally picks 'lm', force 'trf'.
        if robust and (loss != "linear") and method == "lm":
            method = "trf"  # 'lm' supports only linear loss  [oai_citation:4‡docs.scipy.org](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html?utm_source=chatgpt.com)

        # residual function: stack xyz residuals into a 3N vector
        def fun(x):
            ep = ephem0.copy()
            ep[fit_indices] = x
            pred = self.ephem2cart(t_data, ep, compute_velocity=False)
            return (pred - rv_ecef_data[:, :3]).reshape(-1)

        # optional analytic jac is painful here; use finite-difference
        # '2-point' is default; '3-point' is more accurate but more expensive
        from scipy.optimize import least_squares

        res = least_squares(
            fun,
            x0,
            jac="2-point",
            bounds=(lb, ub),
            method=method,
            loss=loss if robust else "linear",
            f_scale=f_scale,
            x_scale=x_scale,
            max_nfev=max_nfev,
            verbose=verbose,
        )  #  [oai_citation:5‡docs.scipy.org](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html?utm_source=chatgpt.com)

        ephem_fit = ephem0.copy()
        ephem_fit[fit_indices] = res.x

        return ephem_fit
