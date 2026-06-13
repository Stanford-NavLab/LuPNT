import numpy as np
import pylupnt as pnt
import matplotlib.pyplot as plt
import copy
from scipy.interpolate import interp1d
import os
from os.path import join, exists
from plotly import graph_objects as go

try:
    from .orbit_utils import rot_moon_ci2pa, convert_ci2pa, convert_pa2ci
    from .orbit_utils import wrapToPi, wrapTo2Pi, rms
except ImportError:
    from orbit_utils import rot_moon_ci2pa, convert_ci2pa, convert_pa2ci
    from orbit_utils import wrapToPi, wrapTo2Pi, rms


class OrbitManager:
    def __init__(self, orbit, dyn, n_period=2, dt=0.1, data_dir=None, overwrite=False):

        self.orbit = orbit
        pnt.set_lupnt_epoch(0.0)

        # time settings
        if orbit == "ELFO" or orbit == "Polar" or orbit == "CLFO" or orbit == "LCRNS":
            sphm = [80, 80]  # Spherical harmonic model degree and order
        else:
            sphm = [8, 8]  # Spherical harmonic model degree and order

        self.n_period = n_period
        self.dt = dt
        self.sphm = sphm

        dyn = self.generate_dynamics()

        print("Setting up orbit: ", orbit)
        self.setup_orbit(orbit)
        self.propagate_orbit(dyn, self.t0_tai, dt, n_period, data_dir=data_dir, overwrite=overwrite)

    def generate_dynamics(self):
        # dynamics
        add_earth = True
        add_sun = True
        sphm = self.sphm

        dyn = pnt.NBodyDynamics()
        dyn.set_integrator(pnt.IntegratorType.RKF45)
        dyn.set_integrator_params(pnt.IntegratorParams(max_iter=20, abstol=1e-12, reltol=1e-12))
        dyn.add_body(pnt.Body.Moon(sphm[0], sphm[1]))
        if add_earth:
            dyn.add_body(pnt.Body.Earth())
        if add_sun:
            dyn.add_body(pnt.Body.Sun())
        dyn.set_frame(pnt.MOON_CI)
        dyn.set_time_step(1.0)  # propagation timestep

        return dyn

    def setup_orbit(self, orbit):

        twobody = True

        if orbit == "ELFO":
            t0_tai = pnt.convert_time(pnt.gregorian_to_time(2025, 1, 1, 12, 0, 0), pnt.UTC, pnt.TAI)

            a = 6541.4e3
            ecc = 0.6
            inc = np.deg2rad(65.5)
            Omega = np.deg2rad(60)
            w = np.deg2rad(90)
            M0 = np.deg2rad(0)

            coe_op = np.array([a, ecc, inc, Omega, w, M0])
            rv0_op = pnt.classical_to_cart(coe_op, pnt.GM_MOON)  # In OP frame
            rv0_ci = pnt.convert_frame(t0_tai, rv0_op, pnt.MOON_OP, pnt.MOON_CI)  # In CI frame

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
            rv0_pa = convert_ci2pa(t0_tai, rv0_ci, rotate_only=True)

        elif orbit == "LNSS":
            t0_tai = pnt.convert_time(pnt.gregorian_to_time(2027, 1, 1, 0, 0, 0), pnt.TDB, pnt.TAI)
            a = 6143.0e3  # meters
            ecc = 0.60
            inc = np.deg2rad(56.2)
            w = np.deg2rad(90.0)
            Omega = np.deg2rad(0.0)
            M0 = np.deg2rad(pnt.true_to_mean_anomaly(np.deg2rad(0.0), ecc))
            coe_op = np.array([a, ecc, inc, Omega, w, M0])
            rv0_op = pnt.classical_to_cart(coe_op, pnt.GM_MOON)
            rv0_ci = pnt.convert_frame(t0_tai, rv0_op, pnt.MOON_OP, pnt.MOON_CI)  # In CI frame
            rv0_pa = convert_ci2pa(t0_tai, rv0_ci, rotate_only=True)

        elif orbit == "LCRNS":
            # https://esc.gsfc.nasa.gov/static-files/LCRNS_Reference_Constellation_White_Paper_03_2025.pdf
            t0_tai = pnt.convert_time(pnt.gregorian_to_time(2027, 3, 1, 0, 0, 0), pnt.UTC, pnt.TAI)

            svid = 0  # 0-4  <--------- chage this to select the satellite

            # a, e, i, Omega, w, M0
            svoe = np.zeros((5, 6))
            svoe[0] = np.array([11315.936501, 0.691982, 59.373229, 321.019197, 92.494031, 0.000000])
            svoe[1] = np.array(
                [11317.948675, 0.691982, 58.951732, 320.997768, 92.505016, 180.000000]
            )
            svoe[2] = np.array(
                [11305.413654, 0.691982, 52.733096, 81.148790, 92.062891, 140.049207]
            )
            svoe[3] = np.array(
                [11326.302154, 0.691982, 52.513419, 81.138818, 92.068945, 195.992393]
            )
            svoe[4] = np.array(
                [11307.882863, 0.691982, 56.310396, 204.889626, 85.444071, 164.007607]
            )
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
            rv0_pa = convert_ci2pa(t0_tai, rv0_ci, rotate_only=True)

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

            rv0_pa = convert_ci2pa(t0_tai, rv0_ci, rotate_only=True)

        elif orbit == "CLFO":
            t0_tai = pnt.convert_time(pnt.gregorian_to_time(2025, 1, 1, 0, 0, 0), pnt.UTC, pnt.TAI)
            a = 6215.0e3
            ecc = 0.0001
            inc = np.deg2rad(39.23)
            Omega = np.deg2rad(0)
            w = np.deg2rad(90)
            M0 = np.deg2rad(60)

            coe_op = np.array([a, ecc, inc, Omega, w, M0])
            rv0_op = pnt.classical_to_cart(coe_op, pnt.GM_MOON)  # In OP frame
            rv0_ci = pnt.convert_frame(t0_tai, rv0_op, pnt.MOON_OP, pnt.MOON_CI)  # In CI frame

        elif orbit == "Halo_S25":
            t0_tai = pnt.convert_time(pnt.gregorian_to_time(2021, 1, 1, 0, 0, 0), pnt.UTC, pnt.TAI)
            rv0_ci = (
                np.array(
                    [
                        -38279.9298,
                        61417.8316,
                        -50632.1007,
                        0.0895177417,
                        0.0538711588,
                        -0.00233218533,
                    ]
                )
                * 1e3
            )
            T_syn = 29.68278536728694 * 86400
            # synodic period
            period = T_syn * 2 / 5
            twobody = False

            coe = pnt.cart_to_classical(rv0_ci, pnt.GM_MOON)
            a = coe[0]
            ecc = coe[1]
            inc = coe[2]

            self.init_phase = 0  # deg

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

            self.init_phase = 0  # deg

        print("Initial position and velocity in CI frame:")
        print(rv0_ci)
        if rv0_ci.shape != (6,):
            raise ValueError("Initial state vector must be 6 elements.")

        if twobody:
            period = 2 * np.pi * np.sqrt(a**3 / pnt.GM_MOON)
            self.init_phase = M0

        print("Semi-major axis [km]:", a / 1000)
        print("Eccentricity:", ecc)
        print("Inclination [deg]:", np.rad2deg(inc))
        print("Period [hr]:", period / 3600)

        coe_pa = pnt.cart_to_classical(rv0_pa, pnt.GM_MOON)
        Omega_deg = np.rad2deg(coe_pa[3])
        f_deg = np.rad2deg(pnt.mean_to_true_anomaly(coe_pa[5], coe_pa[1]))
        if Omega_deg < 0:
            Omega_deg += 360
        if f_deg < 0:
            f_deg += 360

        print(" ")
        print("----------------------------------------")
        print("Initial orbital elements in PA frame:")
        print("a [km]:", coe_pa[0] / 1000)
        print("e:", coe_pa[1])
        print("i [deg]:", np.rad2deg(coe_pa[2]))
        print("Omega [deg]:", Omega_deg)
        print("w [deg]:", np.rad2deg(coe_pa[4]))
        print("f [deg]:", np.rad2deg(pnt.mean_to_true_anomaly(coe_pa[5], coe_pa[1])))
        print("----------------------------------------")
        print(" ")

        # store
        self.t0_tai = t0_tai
        self.rv0_ci = rv0_ci
        self.period = period
        self.twobody = twobody

    def propagate_orbit(self, dyn, t0_tai, dt, n_period, data_dir=None, overwrite=False):

        tlen = n_period * self.period
        n_points = int(tlen / dt) + 1
        tspan = np.linspace(0, tlen, n_points)

        t_tai = t0_tai + tspan
        lent = tspan.size

        # save
        if not exists(data_dir):
            print(f"Creating directory: {data_dir}")
            os.makedirs(data_dir)

        save_file = join(data_dir, f"{self.orbit}_p_{n_period}_dt_{dt}.npz")

        # save to file
        if exists(save_file) and not overwrite:
            print(f"Loading existing data from {save_file}")
            data = np.load(save_file)
            rv_prop_mci = data["rv_ci"]
        else:
            print(f"Propagating orbit and saving to {save_file}")
            # propagate
            rv_prop_mci = dyn.propagate(self.rv0_ci, t_tai)
            # save to file
            np.savez(save_file, rv_ci=rv_prop_mci)

        # convert to orbital elements
        oe_prop = np.zeros(rv_prop_mci.shape)
        for i, rv in enumerate(rv_prop_mci):
            oe_prop[i] = pnt.cart_to_classical(rv, pnt.GM_MOON)

        # Mean anomaly vector
        # if self.twobody:
        #     Mprop = np.rad2deg(oe_prop[:, 5])
        #     Mprop[Mprop <= 0] += 360
        #     M_prop_stack = np.hstack([Mprop, Mprop[1:] + Mprop[-1]])

        #     # fix sudden jumps in M_prop_stack
        #     for i in range(1, M_prop_stack.size):
        #         if M_prop_stack[i] - M_prop_stack[i-1] > 180:
        #             M_prop_stack[i:] -= 360
        #         elif M_prop_stack[i] - M_prop_stack[i-1] < -180:
        #             M_prop_stack[i:] += 360
        # else:
        M_prop_stack = tspan / self.period * 360 + self.init_phase

        # store
        self.rv_prop_mci = rv_prop_mci
        self.oe_prop = oe_prop
        self.M_prop_stack = M_prop_stack
        self.tspan = tspan
        self.t_tai = t_tai
        self.lent = lent

    def plot_coe(self, plot_period=False, plot_M_range=None, fit_poly=False):
        # variables
        lent = self.lent
        tspan = self.tspan
        oe_prop = self.oe_prop

        fig, axs = plt.subplots(4, 2, figsize=(12, 12))

        if plot_M_range is None:
            plot_idx = np.arange(lent)
        else:
            plot_idx = np.where(
                (self.M_prop_stack > plot_M_range[0]) & (self.M_prop_stack < plot_M_range[1])
            )[0]
            print("duration of the plot [min]: ", (tspan[plot_idx[-1]] - tspan[plot_idx[0]]) / 60)

        DEGRAD = 180 / np.pi
        scale = [1, 1, DEGRAD, DEGRAD, DEGRAD, DEGRAD, 1, DEGRAD]
        labels = [
            "a (km)",
            "e",
            "i (deg)",
            "Omega (deg)",
            "w (deg)",
            "Mean Anomaly (deg)",
            "r (km)",
            "u (deg)",
        ]
        for i, ax in enumerate(axs.flatten()):
            if i < 5:
                plot_y = scale[i] * oe_prop[plot_idx, i]
            elif i == 5:
                Mplot = scale[i] * oe_prop[plot_idx, i]
                fplot = np.rad2deg(
                    pnt.mean_to_true_anomaly(np.deg2rad(Mplot), oe_prop[plot_idx, 1])
                )
                Mplot[Mplot <= 0] += 360
                fplot[fplot <= 0] += 360
                for ti in range(1, len(tspan[plot_idx])):
                    if abs(Mplot[ti] - Mplot[ti - 1]) > np.pi:
                        Mplot[ti:] -= 360 * np.sign(Mplot[ti] - Mplot[ti - 1])
                    if abs(fplot[ti] - fplot[ti - 1]) > np.pi:
                        fplot[ti:] -= 360 * np.sign(fplot[ti] - fplot[ti - 1])
                plot_y = Mplot
            elif i == 6:
                plot_y = scale[i] * np.linalg.norm(self.rv_prop_mci[plot_idx, :3], axis=1)
            elif i == 7:
                wplot = scale[4] * oe_prop[plot_idx, 4]
                plot_y = wplot + fplot

            ax.plot(tspan[plot_idx] / 3600, plot_y, "ko-", markersize=1, label=labels[i])

            if fit_poly:
                poly_degs = [1, 2, 3, 4]
                for poly_deg in poly_degs:
                    coeffs = np.polyfit(tspan[plot_idx], plot_y, poly_deg)
                    poly = np.poly1d(coeffs)
                    ax.plot(
                        tspan[plot_idx] / 3600,
                        poly(tspan[plot_idx]),
                        "--",
                        label="Poly fit (degree: {0})".format(poly_deg),
                    )
            ax.set_xlabel("Time (hr)")
            ax.set_ylabel(labels[i])
            ax.grid(True)
            ax.set_xlim([tspan[plot_idx[0]] / 3600, tspan[plot_idx[-1]] / 3600])
            ax.legend()
            ax.set_title(
                [
                    "Semi-major axis",
                    "Eccentricity",
                    "Inclination",
                    "RAAN",
                    "Argument of perigee",
                    "Mean Anomaly",
                    "Radial Distance",
                    "Argument of Latitude",
                ][i]
            )

        if plot_period:
            a = np.mean(oe_prop[:, 0])
            tend = tspan[-1]
            t_orb = tspan[0]
            while t_orb <= tend:
                for ax in axs.flatten():
                    ax.axvline(
                        x=t_orb / 3600,
                        color="black",
                        linestyle="--",
                        linewidth=1.5,
                    )
                t_orb += self.period

        plt.tight_layout()

    def plot_orbit(self, inv=600, figname=None, use_black_back=False):
        # 3d plot of the orbit
        t_tai = self.t_tai[::inv]
        rv_prop_mci = self.rv_prop_mci[::inv, :]
        orbit = self.orbit

        # convert to moon fixed frame
        rv_prop_bf = convert_ci2pa(t_tai, rv_prop_mci, rotate_only=True)

        fig = go.Figure()

        if orbit == "DRO":
            rv_prop_eci = rv_prop_mci + pnt.get_body_pos_vel(
                t_tai, pnt.EARTH, pnt.MOON, pnt.Frame.ECI
            )  # (S - M) + (M - E) = S - E
            lunar_orbit_eci = pnt.get_body_pos_vel(t_tai, pnt.MOON, pnt.Frame.ECI)

            print(lunar_orbit_eci.shape)

            lunar_orbit_pa = convert_ci2pa(t_tai, lunar_orbit_eci, rotate_only=True)
            rv_prop_pa = convert_ci2pa(t_tai, rv_prop_eci, rotate_only=True)

            orbits = np.zeros((2, t_tai.size, 6))
            orbits[0] = lunar_orbit_eci
            orbits[1] = rv_prop_eci
            pnt.plot.plot_orbits(fig, orbits, color=["gray", "blue"])

        elif orbit == "NRHO" or orbit == "Halo_S25":
            pnt.plot.plot_orbits(fig, rv_prop_bf, color="blue")
        else:
            pnt.plot.plot_orbits(fig, rv_prop_mci, color="blue")

        pnt.plot.plot_body(
            fig,
            pnt.MOON,
            size_factor=2,
            alpha=0.5,
        )

        if orbit == "LCRNS":
            zoom = 3.0
            az = -45
            el = 20
        elif orbit == "Polar":
            az = -45
            el = 20
            zoom = 3.0
        elif orbit == "NRHO":
            az = -45
            el = 30
            zoom = 3.5
        else:
            az = -45
            el = 20
            zoom = 2.5

        pnt.plot.set_view(fig, az, el, zoom)
        fig.update_layout(showlegend=True, width=400, height=400)

        if use_black_back:
            # make background black
            fig.update_layout(paper_bgcolor="black", plot_bgcolor="black")
            # delete axis
            fig.update_layout(
                scene=dict(
                    xaxis=dict(
                        showbackground=False, visible=False, showticklabels=False, showgrid=False
                    ),
                    yaxis=dict(
                        showbackground=False, visible=False, showticklabels=False, showgrid=False
                    ),
                    zaxis=dict(
                        showbackground=False, visible=False, showticklabels=False, showgrid=False
                    ),
                )
            )
            # delete labels
            fig.update_layout(scene=dict(xaxis_title="", yaxis_title="", zaxis_title=""))

        if figname is not None:
            fig.write_image(figname)
        fig.show()

    def extract_orbit_prop(
        self,
        center_M,
        fit_min,
        n_points,
        use_chebyshev_sample=False,
        plot_points=False,
        p_cheby=1.0,
    ):
        # variables
        fit_t = fit_min * 60 / 2
        dM = fit_t / self.period * 360
        start_M = center_M - dM  # start of the orbit (mean anomaly)
        end_M = center_M + dM  # end of the orbit (mean anomaly)
        dM_start = start_M - self.M_prop_stack[0]
        dM_end = end_M - self.M_prop_stack[0]
        if dM_start < 0:
            start_M += 360
            dM_start += 360
        if dM_end < dM_start:
            end_M += 360
            dM_end += 360

        if use_chebyshev_sample:  # sample by Chebyshev points
            chebyshev_points = np.cos(np.pi * (np.arange(0, n_points + 1) / n_points))[::-1]
            chebyshev_points = np.sign(chebyshev_points) * np.power(
                np.abs(chebyshev_points), p_cheby
            )
        else:  # sample by mean anomaly
            M_points = np.linspace(start_M, end_M, n_points + 1)

        if plot_points:
            print("M_points: ", M_points)

        t_data = []
        dM_prev = 0.0
        for mi, M in enumerate(M_points):
            dM_point = M - self.M_prop_stack[0]
            if dM_point < 0:
                dM_point += 360
            if mi >= 1:
                if dM_point - dM_prev > 180:
                    dM_point -= 360
                elif dM_point - dM_prev < -180:
                    dM_point += 360
            t_data.append(dM_point * self.period / 360)
            dM_prev = copy.deepcopy(dM_point)

        # propagate orbit
        t_data = np.array(t_data)
        if plot_points:
            print("t_data (hr): ", t_data / 3600)
        tai_data = self.t0_tai + t_data

        dyn = self.generate_dynamics()
        rvi = dyn.propagate(self.rv_prop_mci[0], tai_data)

        rvbf = convert_ci2pa(tai_data, rvi, rotate_only=True)
        rvbf_w = convert_ci2pa(tai_data, rvi, rotate_only=False)

        if plot_points:
            fig, axes = plt.subplots(3, 2, figsize=(8, 8))
            coe_b = np.zeros_like(rvbf)
            f_prev = -np.inf
            for k in range(rvbf.shape[0]):
                coe = pnt.cart_to_classical(rvbf[k], pnt.GM_MOON)
                coe[5] = pnt.mean_to_true_anomaly(coe[5], coe[1])
                coe[2:] = np.rad2deg(coe[2:])
                if coe[5] < 0:
                    coe[5] += 360
                if coe[5] < f_prev:
                    coe[5] += 360
                coe_b[k] = coe
                f_prev = copy.deepcopy(coe[5])

            for i, ax in enumerate(axes.flatten()):
                ax.plot(
                    (tai_data - tai_data[0]) / 60, coe_b[:, i], "bo-", markersize=1, linewidth=1
                )
                # ax.plot((t_data - t_data[0])/60, rvbf[:, i], 'bo-')
                ax.set_xlabel("Time (min)")
                ax.set_ylabel(["a", "e", "i [deg]", "lambda [deg]", "w [deg]", "f [deg]"][i])
                ax.grid(True)

            plt.suptitle(
                "Orbit from M = {0:.1f} deg to {1:.1f} deg ({2} min)".format(
                    start_M, end_M, fit_min
                )
            )
            plt.tight_layout()

        # bring to 0
        t_data -= t_data[0]

        return t_data, rvi, rvbf, rvbf_w

    def extract_orbit_interp(
        self,
        center_M,
        fit_min,
        n_points,
        use_chebyshev_sample=False,
        plot_points=False,
        p_cheby=1.0,
        debug=False,
    ):
        # variables
        fit_t = fit_min * 60 / 2
        dM = fit_t / self.period * 360
        start_M = center_M - dM  # start of the orbit (mean anomaly)
        end_M = center_M + dM  # end of the orbit (mean anomaly)
        dM_start = start_M - self.M_prop_stack[0]
        dM_end = end_M - self.M_prop_stack[0]
        if dM_start < 0:
            start_M += 360
            dM_start += 360
        if dM_end < dM_start:
            end_M += 360
            dM_end += 360

        # find the closest index
        margin = 600
        start_idx = max(np.argmin(np.abs(self.M_prop_stack - start_M)) - margin, 0)
        end_idx = min(np.argmin(np.abs(self.M_prop_stack - end_M)) + margin, self.lent)
        # convert true anomaly to time
        tai_data_span = self.t_tai[start_idx:end_idx]
        t_data_span = self.tspan[start_idx:end_idx]
        rvi_span = self.rv_prop_mci[start_idx:end_idx]

        rvbf_span = convert_ci2pa(tai_data_span, rvi_span, rotate_only=True)
        rvbf_w_span = convert_ci2pa(tai_data_span, rvi_span, rotate_only=False)

        # evaluation points
        t_mid = (t_data_span[0] + t_data_span[-1]) / 2
        t_start = t_mid - fit_t
        t_end = t_mid + fit_t
        t_data = np.linspace(t_start, t_end, n_points + 1)
        if use_chebyshev_sample:  # sample by Chebyshev points
            chebyshev_points = np.cos(np.pi * (np.arange(0, n_points + 1) / n_points))[::-1]
            chebyshev_points = np.sign(chebyshev_points) * np.power(
                np.abs(chebyshev_points), p_cheby
            )
            t_data = 0.5 * (t_end - t_start) * chebyshev_points + 0.5 * (t_start + t_end)

        # interpolate
        tai_data = interp1d(t_data_span, tai_data_span, kind="linear", fill_value="extrapolate")(
            t_data
        )
        rvi = interp1d(t_data_span, rvi_span, axis=0, kind="cubic", fill_value="extrapolate")(
            t_data
        )
        rvbf = interp1d(t_data_span, rvbf_span, axis=0, kind="cubic", fill_value="extrapolate")(
            t_data
        )
        rvbf_w = interp1d(t_data_span, rvbf_w_span, axis=0, kind="cubic", fill_value="extrapolate")(
            t_data
        )

        np.set_printoptions(precision=5, suppress=True)

        if debug:
            print("start_idx (M={0:.1f}) : {1}".format(start_M, start_idx))
            print("end_idx   (M={0:.1f}) : {1}".format(end_M, end_idx))
            print("querty time (hr)   : ", t_data_span / 3600)
            print("dt          (s)    : ", t_data_span[1] - t_data_span[0])
            print("query tspan (min)  : ", (t_data_span[-1] - t_data_span[0]) / 60)
            print("sampled time (hr)  : ", t_data / 3600)
            print("sampled tspan (min): ", (t_data[-1] - t_data[0]) / 60)
            print("number of points   : ", n_points)
            print(" ")

        if plot_points:
            fig, axes = plt.subplots(3, 2, figsize=(8, 8))
            coe_b = np.zeros_like(rvbf)
            for k in range(rvbf.shape[0]):
                coe_b[k] = pnt.cart_to_classical(rvbf[k], pnt.GM_MOON)
                coe_b[k, 2:] = np.rad2deg(coe_b[k, 2:])

            for i, ax in enumerate(axes.flatten()):
                ax.plot(
                    (tai_data - tai_data[0]) / 60, coe_b[:, i], "bo-", markersize=1, linewidth=1
                )
                # ax.plot((t_data - t_data[0])/60, rvbf[:, i], 'bo-')
                ax.set_xlabel("Time (min)")
                ax.set_ylabel(
                    ["a [km]", "e []", "i [deg]", "lambda [deg]", "w [deg]", "M [deg]"][i]
                )
                ax.grid(True)

            plt.suptitle(
                "Orbit from M = {0:.1f} deg to {1:.1f} deg ({2} min)".format(
                    start_M, end_M, fit_min
                )
            )
            plt.tight_layout()

        # bring to 0
        # t_data -= t_data[0]

        return t_data, rvi, rvbf, rvbf_w

    def extract_orbit(
        self,
        center_M,
        fit_min,
        method="prop",
        n_points=None,
        use_chebyshev_sample=False,
        plot_points=False,
        p_cheby=1.0,
        debug=False,
    ):
        # variables
        if method == "prop":
            t_data, rvi, rvbf, rvbf_w = self.extract_orbit_prop(
                center_M, fit_min, n_points, use_chebyshev_sample, plot_points, p_cheby=p_cheby
            )
        elif method == "interp":
            t_data, rvi, rvbf, rvbf_w = self.extract_orbit_interp(
                center_M,
                fit_min,
                n_points,
                use_chebyshev_sample,
                plot_points,
                p_cheby=p_cheby,
                debug=debug,
            )
        else:
            raise ValueError("Method not supported: ", method)

        return t_data, rvi, rvbf, rvbf_w
