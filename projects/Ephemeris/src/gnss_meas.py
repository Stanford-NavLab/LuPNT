import numpy as np
import pandas as pd
import pylupnt as pnt
import plotly.graph_objects as go

try:
    from .odmeas import ODMeas
except ImportError:
    from odmeas import ODMeas


class GNSSReceiverParam:
    def __init__(self):
        # GNSS receiver parameters
        self.Tsys = 190.0  # System noise temp [K]
        self.Ae = 0.0  # Attenuation due to atmosphere (should be negative) [dB]
        self.Nf = -2.85  # Noise figure of receiver/LNA [dB]
        self.L = -0.16  # Receiver implementation, A/D conversion losses [dB]
        self.As = 0.0  # System losses, in front of LNA [dB]
        self.CN0_threshold = 15.0  # CN0 threshold [dB-Hz]

        # GNSS Receiver chip parameter
        self.Bp = 5.0  # Carrier loop noise bandwidth [Hz]
        self.T = 20e-3  # Tracking loop integration time [s]
        self.b = 2.0  # normalized bandwidth [Hz]
        self.Bn = 0.2  # Code loop noise bandwidth [Hz]
        self.D = 0.3  # Early-to-late correlator spacing (chips)


gnss_freq_map = {
    "L1": 1575.42e6,
    "L2": 1227.60e6,
    "L5": 1176.45e6,
    "E1": 1575.42e6,
    "E6": 1278.75e6,
    "E5": 1191.795e6,
    "E5a": 1176.45e6,
    "E5b": 1207.14e6,
}

gnss_rc_map = {
    "L1": 1.023e6,
    "L2": 0.5115e6,
    "L5": 10.23e6,
    "E1": 1.023e6,
    "E6": 0.5115e6,
    "E5": 10.23e6,
    "E5a": 10.23e6,
    "E5b": 10.23e6,
}


class GNSSMeas(ODMeas):
    def __init__(self, t_tai, rv_m2sc_ci, gnss_consts=["GPS", "GALILEO"]):
        super().__init__("GPS")
        self.n_meas = 0

        # dynamics with J2
        dyn_gnss = pnt.NBodyDynamics(pnt.IntegratorType.RKF45)
        dyn_gnss.set_integrator_params(
            pnt.IntegratorParams(max_iter=20, abstol=1e-10, reltol=1e-10)
        )
        dyn_gnss.add_body(pnt.Body.Earth(2, 0))  # Earth
        dyn_gnss.set_frame(pnt.ECI)
        dyn_gnss.set_time_step(600)
        dyn_gnss.set_print_progress(False)
        self.dyn_gnss = dyn_gnss

        # initialize the gnss data
        self.t_tai = t_tai
        self.rv_m2sc_ci = rv_m2sc_ci  # spacecraft position in Moon CI frame

        # GNSS receiver parameters
        self.gnssr_param = GNSSReceiverParam()

        self.setup_gnss(gnss_consts)

    def compute_cn0(
        self,
        t_tai,
        rv_m2sc_ci,
        N_gnss,
        rv_gnss_ci,
        ex_gnss,
        ey_gnss,
        ez_gnss,
        gnss_antennas,
        prn_gnss,
        P_tx,
    ):

        # numbers
        N_sc = rv_m2sc_ci.shape[0]  # number of spacecraft
        N_t = rv_m2sc_ci.shape[1]  # number of time steps

        # Visibility
        vis_sc2sc = np.ones((N_sc, N_sc, N_t), dtype=bool)
        vis_sc2gps = np.ones((N_sc, N_gnss, N_t), dtype=bool)
        dist_sc2sc = np.zeros((N_sc, N_sc, N_t))
        dist_sc2gps = np.zeros((N_sc, N_gnss, N_t))
        phi_gps2sc = np.zeros((N_sc, N_gnss, N_t))
        phi_sc2gps = np.zeros((N_sc, N_gnss, N_t))
        theta_gps2sc = np.zeros((N_sc, N_gnss, N_t))
        G_tx = np.zeros((N_sc, N_gnss, N_t))
        G_rx = np.zeros((N_sc, N_gnss, N_t))

        # sun direction
        # rv_e2s_eci = pnt.get_body_pos_vel(t_tai, pnt.EARTH, pnt.SUN, pnt.ECI)
        rv_m2e_ci = pnt.get_body_pos_vel(t_tai, pnt.MOON, pnt.EARTH, pnt.MOON_CI)

        for i in range(N_sc):
            for j in range(i + 1, N_sc):
                vis_sc2sc[i, j] &= self.compute_visibility(
                    rv_m2sc_ci[i, :, :3], rv_m2sc_ci[j, :, :3], pnt.R_MOON
                )
                vis_sc2sc[j, i] &= vis_sc2sc[i, j]
                dist_sc2sc[i, j] = np.linalg.norm(
                    rv_m2sc_ci[i, :, :3] - rv_m2sc_ci[j, :, :3], axis=1
                )
                dist_sc2sc[j, i] = dist_sc2sc[i, j]

            for j in range(N_gnss):
                # moon blockage
                vis_sc2gps[i, j] &= self.compute_visibility(
                    rv_m2sc_ci[i, :, :3], rv_gnss_ci[j, :, :3], pnt.R_MOON
                )

                # earth blockage
                vis_sc2gps[i, j] &= self.compute_visibility(
                    rv_m2sc_ci[i, :, :3],
                    rv_gnss_ci[j, :, :3],
                    pnt.R_EARTH,
                    rv_m2e_ci[:, :3],
                )

                dist_sc2gps[i, j] = np.linalg.norm(
                    rv_m2sc_ci[i, :, :3] - rv_gnss_ci[j, :, :3], axis=1
                )

                u_gps2sc = pnt.normalize(rv_m2sc_ci[i, :, :3] - rv_gnss_ci[j, :, :3])
                phi_gps2sc[i, j] = np.arccos(np.clip(np.sum(u_gps2sc * ez_gnss[j], axis=-1), -1, 1))
                theta_gps2sc[i, j] = np.arctan2(
                    np.sum(u_gps2sc * ey_gnss[j], axis=-1),
                    np.sum(u_gps2sc * ex_gnss[j], axis=-1),
                )

                u_sc2gps = -u_gps2sc
                u_sc2e = pnt.normalize(rv_m2e_ci[:, :3] - rv_m2sc_ci[i, :, :3])
                phi_sc2gps[i, j] = np.arccos(np.clip(np.sum(u_sc2e * u_sc2gps, axis=-1), -1, 1))

                G_tx[i, j] = gnss_antennas[prn_gnss[j]].compute_gain(
                    theta_gps2sc[i, j], phi_gps2sc[i, j]
                )
                G_rx[i, j] = self.rx_antenna.compute_gain(0, phi_sc2gps[i, j])

                G_tx[i, j][~vis_sc2gps[i, j]] = np.nan
                G_rx[i, j][~vis_sc2gps[i, j]] = np.nan

            vis_sc2sc[i, i, :] = False

        # Link budget
        freq = 1575.42e6  # [Hz] GPS L1 frequency
        L_ad = 0.0  # [dB] A/D converter loss
        L_atm = 0.0  # [dB] Atmospheric loss
        Nf = 2  # [dB] Noise figure
        Tsys = 113  # [K] System noise temperature
        # [dB] Free space loss
        L_fs = 20 * np.log10((4 * np.pi * dist_sc2gps) / (pnt.C / freq))
        # [dB] Transmitter antenna gain
        kb_db = 10 * np.log10(1.38e-23)  # [dB] Boltzmann constant

        CN0 = P_tx + G_tx + G_rx - L_atm - L_fs - L_ad - Nf - kb_db - 10 * np.log10(Tsys)

        return CN0

    def compute_visibility(
        self, r1: np.ndarray, r2: np.ndarray, R_body: float, r_body: np.ndarray = None
    ) -> np.ndarray:
        if r_body is None:
            r_body = np.zeros(3)
        r = r2 - r1
        r_norm = np.linalg.norm(r, axis=1)
        r1body = r1 - r_body
        r1body_norm = np.linalg.norm(r1body, axis=1)
        dot = np.einsum("ij,ij->i", -r1body, r)
        theta1 = np.arccos(np.clip(dot / r_norm / r1body_norm, -1, 1))
        theta2 = np.arcsin(np.clip(R_body / r1body_norm, -1, 1))
        visibility = np.ones(len(r1), dtype=bool)
        visibility[(theta1 < theta2) & (r_norm > r1body_norm)] = False
        return visibility

    def setup_gnss(self, gnss_consts=["GPS"]):
        self.gnss_const = gnss_const
        self.n_meas = 0
        t_tai = self.t_tai

        N_t = t_tai.shape[0]
        self.gnss_consts = gnss_consts

        tle_files = {
            "GPS": "gps_2025_01_01",
            "GALILEO": "galileo_2025_01_01",
            "QZSS": "qzss_2025_01_01",
        }

        # sun direction
        rv_e2s_eci = pnt.get_body_pos_vel(t_tai, pnt.EARTH, pnt.SUN, pnt.ECI)

        # variables to store the results
        self.prn_gnss = {}
        self.rv_gnss_eci = {}
        self.rv_gnss_ci = {}
        self.true_range = {}
        self.true_rangerate = {}
        self.true_carrier_phase = {}
        self.sigma_range = {}
        self.sigma_rangerate = {}
        self.sigma_carrier_phase = {}
        self.cn0 = {}
        self.vis_gnss = {}
        self.N_gnss = {}

        for gnss_const in gnss_consts:
            tles = pnt.TLE.from_file(tle_files[gnss_const])
            N_gnss = len(tles)
            rv_gnss_eci = np.zeros((N_gnss, N_t, 6))
            rv_gnss_ci = np.zeros((N_gnss, N_t, 6))

            prns = np.zeros(N_gnss, dtype=int)
            for i in range(N_gnss):
                print(f"Propagating {gnss_const} satellite {i+1}/{N_gnss}")
                coe0_gps = pnt.tle2classical(tles[i], pnt.GM_EARTH)
                rv0_gps_eci = pnt.classical_to_cart(coe0_gps, pnt.GM_EARTH)
                rv_gnss_eci[i] = self.dyn_gnss.propagate(rv0_gps_eci, tles[i].epoch_tai, t_tai)
                rv_gnss_ci[i] = pnt.convert_frame(t_tai, rv_gnss_eci[i], pnt.ECI, pnt.MOON_CI)
                prns[i] = tles[i].prn

            # attitde computation
            e_gnss2e = pnt.normalize(-rv_gnss_eci[:, :, :3])
            e_gnss2s = pnt.normalize(rv_e2s_eci[None, :, :3] - rv_gnss_eci[:, :, :3])
            ez_gnss = e_gnss2e
            ey_gnss = pnt.cross_norm(e_gnss2e, e_gnss2s)
            ex_gnss = pnt.cross_norm(ey_gnss, ez_gnss)

            # antennas
            if gnss_const == "GPS":
                df = pd.read_csv(pnt.find_file("gps_table.csv"))
                gps_antenna_names = df.set_index("PRN")["LM_File"].to_dict()
                # for gps that LM file is not available, use the ACE_file instead
                for prn in gps_antenna_names.keys():
                    if type(gps_antenna_names[prn]) == float:
                        gps_antenna_names[prn] = df.loc[prn - 1, "ACE_File"]
                print(gps_antenna_names)
                antennas = {k: pnt.Antenna(v) for k, v in gps_antenna_names.items()}
            elif gnss_const == "GALILEO":
                antennas = {k: pnt.Antenna("Galileo_E1") for k in prns}
            elif gnss_const == "QZSS":
                qzss_names = ["1R", "02", "03", "04", "05", "06", "07"]
                antennas = {
                    k: pnt.Antenna("QZSS_" + qzss_names[k - 1] + "_L1") for k in range(1, 5)
                }

            # Compute the CN0
            CN0 = self.compute_cn0(
                self.rv_m2sc_ci,
                N_gnss,
                rv_gnss_ci,
                ex_gnss,
                ey_gnss,
                ez_gnss,
                antennas,
                prns,
                pnt.P_TX_GNSS,
            )
            vis_gnss = CN0 > pnt.CN0_THRESHOLD

            # compute the true range, range rate, and carrier phase
            sigma_range = np.zeros((N_gnss, N_t))
            sigma_rangerate = np.zeros((N_gnss, N_t))
            sigma_carrier_phase = np.zeros((N_gnss, N_t))
            sigma_range[vis_gnss] = self.compute_gnss_pseudorange_noise(
                CN0[vis_gnss], self.gnssr_param, "L1"
            )
            sigma_rangerate[vis_gnss] = self.compute_gnss_pseudorangerate_noise(CN0[vis_gnss], "L1")
            sigma_carrier_phase[vis_gnss] = self.compute_gnss_carrier_phase_noise(
                CN0[vis_gnss], "L1"
            )

            # storage
            self.prn_gnss[gnss_const] = prns
            self.rv_gnss_eci[gnss_const] = rv_gnss_eci
            self.rv_gnss_ci[gnss_const] = rv_gnss_ci
            self.cn0[gnss_const] = CN0
            self.vis_gnss[gnss_const] = vis_gnss

    def compute_gnss_pseudorange_noise(CN0_dB, gnssr_param, freq):
        CN0 = 10 ** (CN0_dB / 10)
        Bn = gnssr_param.Bn
        Rc = gnss_rc_map(freq)  # Chip rate in Hz
        Bfe = gnssr_param.b * Rc
        T = gnssr_param.T
        D = gnssr_param.D
        Tc = 1.0 / Rc
        C = 299792.458  # Speed of light in km/s

        sigma = np.zeros_like(CN0)

        case1 = D >= (np.pi * Rc / Bfe)
        case2 = (D > (Rc / Bfe)) & (~case1)
        case3 = ~case1 & ~case2

        if case1:
            sigma = np.sqrt(Bn / (2.0 * CN0) * D * (1.0 + 2.0 / (T * CN0 * (2 - D))))
        if case2:
            tmp1 = Bn / (2.0 * CN0)
            tmp2 = 1.0 / (Bfe * Tc) + Bfe * Tc / (np.pi - 1) * (D - 1.0 / (Bfe * Tc)) ** 2
            tmp3 = 1.0 + 2.0 / (T * CN0 * (2 - D))
            sigma = np.sqrt(tmp1 * tmp2 * tmp3)
        if case3:
            sigma = np.sqrt(Bn / (2.0 * CN0) * (1.0 / (Bfe * Tc)) * (1.0 + 1.0 / (T * CN0)))

        return sigma * (C * Tc)

    def compute_gnss_pseudorangerate_noise(self, CN0_dB, freq):
        gnssr_param = self.gnssr_param
        C = 299792.458  # Speed of light in km/s
        f = gnss_freq_map(freq)
        lambda_ = C / f  # Wavelength of the GNSS signal

        CN0 = 10 ** (CN0_dB / 10)
        F = 2
        Bn = gnssr_param.Bn
        T = gnssr_param.T

        return lambda_ / (2 * np.pi * T) * np.sqrt(4 * F * Bn / CN0 * (1 + 1.0 / (T * CN0)))

    def compute_gnss_carrier_phase_noise(self, CN0_dB, freq):
        CN0 = 10 ** (CN0_dB / 10)
        C = 299792.458  # Speed of light in km/s
        f = gnss_freq_map(freq)
        lambda_ = C / f  # Wavelength of the GNSS signal

        gnssr_param = self.gnssr_param
        Bp = gnssr_param.Bp
        T = gnssr_param.T

        return lambda_ / (2 * np.pi) * np.sqrt(Bp / CN0 * (1.0 + 1.0 / (2 * T * CN0)))

    def plot_gnss_orbit(self, plot_orientation=True):
        tickvals = np.arange(-30, 31, 10) * 1e3
        fig = go.Figure()

        pnt.plot.plot_body(fig, pnt.EARTH, size_factor=5)

        # gps
        for gnss_const in self.gnss_consts:
            pnt.plot.plot_orbits(
                fig, self.rv_gps_eci[:, : int(pnt.SECS_DAY / dt) : 10], t=t, color="lightgray"
            )

            if plot_orientation:
                for i in range(N_gps):
                    pnt.plot.plot_frame(
                        fig,
                        rv_gps_eci[i, t, :3],
                        np.vstack((ex_gps[i, t], ey_gps[i, t], ez_gps[i, t])),
                        length=pnt.R_EARTH,
                        width=5,
                        tip=5,
                    )
                pnt.plot.plot_arrow3(
                    fig,
                    np.zeros(3),
                    pnt.normalize(rv_e2s_eci[t]),
                    length=3 * pnt.R_EARTH,
                    width=5,
                    color="orange",
                    tip=10,
                )
