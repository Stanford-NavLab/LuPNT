from re import match
import numpy as np
import pandas as pd
import pylupnt as pnt
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from datetime import datetime, timedelta
from tqdm import tqdm
import os
import requests
from datetime import datetime
import copy
from src.interfaces.gnss_file_loader import SP3Loader, BRDCLoader
from src.interfaces.gnss_utils import datetime_to_tai
from pylupnt import Logger


def normalize(v):
    norm = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / norm


def cross_norm(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Compute the norm of the cross product of two vectors

    Args:
        a (np.ndarray): first vector
        b (np.ndarray): second vector
    Returns:
        np.ndarray: norm of the cross product
    """
    cross = np.cross(a, b, axis=-1)
    return cross / np.linalg.norm(cross, axis=-1)[..., np.newaxis]


def linkbudget_receiver(freq, dist_sc2gps, P_tx, G_tx, G_rx):
    # NaviMoon Receiver
    # https://navisp.esa.int/uploads/files/documents/NaviMIOD_Final_Presentation.pdf
    # LuGRE receiver
    # https://ntrs.nasa.gov/api/citations/20220010106/downloads/AAS_2022_LuGRE_Analysis_STRIVES_v2.pdf
    freq_Hz = gnss_freq_map[freq]
    L_ad = 0.6  # [dB] A/D converter loss
    L_pol = 1.0  # [dB] Polarization loss
    L_atm = 0.0  # [dB] Atmospheric loss
    Teff = 167.98  # [K] Effective noise temperature
    # [dB] Free space loss
    L_fs = 20 * np.log10((4 * np.pi * dist_sc2gps) / (pnt.C / freq_Hz))
    # [dB] Transmitter antenna gain
    kb_db = 10 * np.log10(1.38e-23)  # [dB] Boltzmann constant

    CN0 = P_tx + G_tx + G_rx - L_atm - L_fs - L_ad - L_pol - kb_db - 10 * np.log10(Teff)

    return CN0


def block_power(block_name, freq):
    # https://ntrs.nasa.gov/api/citations/20220010106/downloads/AAS_2022_LuGRE_Analysis_STRIVES_v2.pdf
    if freq == "L1":
        freq_offset = 0.0
    elif freq == "L5":
        # L5 frequency has a higher transmit power
        # https://insidegnss.com/gps-l5-signal-goes-on-the-air-april-10/
        if block_name in ["IIF", "III"]:
            freq_offset = 3.0
        elif block_name in ["IIR", "IIR_M"]:
            freq_offset = -100.0  # IIR and IIR-M does not transmit L5 signal

    if block_name in ["IIR"]:
        return 17.3 + freq_offset
    elif block_name in ["IIR_M"]:
        return 18.8 + freq_offset
    elif block_name in ["IIF"]:
        return 16.2 + freq_offset
    elif block_name in ["III"]:
        return 18.8 + freq_offset
    else:
        return 14.0  # default


def get_fault_prns(gnss_const, use_faults=False):
    # Return a list of faulty PRNs for a given GNSS constellation
    if not use_faults:
        return []

    if gnss_const == "GPS":
        # https://www.navcen.uscg.gov/gps-constellation
        # Got information for 3/4/2025 using wayback machine
        # https://web.archive.org/web/20250304172707/https://www.navcen.uscg.gov/gps-constellation
        return [6, 8, 21]
    elif gnss_const == "GALILEO":
        # https://www.gsc-europa.eu/system-service-status/constellation-information
        # got information for 2025/3/10 using wayback machine
        # https://web.archive.org/web/20250310172707/https://www.gsc-europa.eu/system-service-status/constellation-information
        return [1, 14, 18, 20, 22]
    elif gnss_const == "QZSS":
        return []
    else:
        return []


class GNSSReceiverParam:
    def __init__(self):
        # GNSS Receiver chip parameter
        # https://www.mdpi.com/1424-8220/16/3/347?utm_source=chatgpt.com
        self.Bp = 1.0  # Carrier loop noise bandwidth [Hz]  (higher Bp (5.0 Hz) for closer distance to Earth)
        self.T = 20e-3  # Tracking loop integration time [s]
        self.b = 2.0  # normalized bandwidth [Hz]
        self.Bn = 0.7  # Code loop noise bandwidth [Hz]
        self.Bf = 0.2  # Frequency loop noise bandwidth [Hz]
        self.D = 0.1  # Early-to-late correlator spacing (chips, correlator spacing should be larger than 0.25 chips)


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


def get_gps_table():
    data = [
        [1, 63, "IIF", "", "SVN63_ACE"],
        [2, 61, "IIR", "SVN61_LM", "SVN61_ACE"],
        [3, 69, "IIF", "", "SVN69_ACE"],
        [4, 74, "III", "SVN74_LM", ""],
        [5, 50, "IIR_M", "SVN50_LM", "SVN50_ACE"],
        [6, 67, "IIF", "", "SVN67_ACE"],
        [7, 48, "IIR_M", "SVN48_LM", "SVN48_ACE"],
        [8, 72, "IIF", "", "SVN72_ACE"],
        [9, 68, "IIF", "", "SVN68_ACE"],
        [10, 73, "IIF", "", "SVN73_ACE"],
        [11, 78, "III", "SVN78_LM", ""],
        [12, 58, "IIR_M", "SVN58_LM", "SVN58_ACE"],
        [13, 43, "IIR", "SVN43_LM", "SVN43_ACE"],
        [14, 77, "III", "SVN77_LM", ""],
        [15, 55, "IIR_M", "SVN55_LM", "SVN55_ACE"],
        [16, 56, "IIR", "SVN56_LM", "SVN56_ACE"],
        [17, 53, "IIR_M", "SVN53_LM", "SVN53_ACE"],
        [18, 75, "III", "SVN75_LM", ""],
        [19, 59, "IIR", "SVN59_LM", "SVN59_ACE"],
        [20, 51, "IIR", "SVN51_LM", "SVN51_ACE"],
        [21, 45, "IIR", "SVN45_LM", "SVN45_ACE"],
        [22, 47, "IIR", "SVN47_LM", "SVN47_ACE"],
        [23, 76, "III", "SVN76_LM", ""],
        [24, 65, "IIF", "", "SVN65_ACE"],
        [25, 62, "IIF", "", "SVN62_ACE"],
        [26, 71, "IIF", "", "SVN71_ACE"],
        [27, 66, "IIF", "", "SVN66_ACE"],
        [28, 79, "III", "SVN78_LM", ""],  # Note: LM file is SVN78_LM
        [29, 57, "IIR_M", "SVN57_LM", "SVN57_ACE"],
        [30, 64, "IIF", "", "SVN64_ACE"],
        [31, 52, "IIR_M", "SVN52_LM", "SVN52_ACE"],
        [32, 70, "IIF", "", "SVN70_ACE"],
    ]
    df = pd.DataFrame(data, columns=["PRN", "SVN", "blockName", "LM_File", "ACE_File"])
    return df


class GNSSMeas:
    def __init__(self, t_tai, rv_m2sc_ci, basepath, consider_faults=False, rv_e2sc_ecef=None):
        self.n_meas = 0
        self.consider_faults = consider_faults

        # initialize the gnss data
        self.t_tai = t_tai
        self.tspan = t_tai - t_tai[0]  # time span in seconds
        self.rv_m2sc_ci = rv_m2sc_ci  # spacecraft position in Moon CI frame
        self.basepath = basepath
        self.N_t = t_tai.shape[0]  # number of time steps
        self.N_s = rv_m2sc_ci.shape[0]  # number of spacecraft

        # convert to gregorian time
        # str_t = pnt.time2gregorian_string(t_tai[0])
        str_t = pnt.time_to_gregorian_string(t_tai[0])
        # yyyy-mm-ddThh:mm:ss.ssssss
        # parse by T
        split_strT = str_t.split("T")[0].split("-")
        self.year = int(split_strT[0])
        self.month = int(split_strT[1])
        self.day = int(split_strT[2])

        self.gps_datetime = None  # GPS datetime for SP3 file loading

        # _, self.df_gnss = get_gnss_table(datetime(self.year, self.month, self.day))

        # convert to ecef
        if rv_e2sc_ecef is not None:
            self.rv_e2sc_ecef = rv_e2sc_ecef
        else:
            self.rv_e2sc_ecef = np.zeros_like(rv_m2sc_ci)
            N_sc = rv_m2sc_ci.shape[0]  # number of spacecraft
            t_tdb = pnt.convert_time(t_tai, pnt.TAI, pnt.TDB)
            for i in range(N_sc):
                self.rv_e2sc_ecef[i] = pnt.convert_frame(
                    t_tdb, rv_m2sc_ci[i], pnt.MOON_CI, pnt.ECEF
                )

        # dyn_prop = pnt.CartesianTwoBodyDynamics(pnt.GM_EARTH, pnt.RKF45)
        dyn_prop = pnt.CartesianTwoBodyDynamics(pnt.GM_EARTH)
        dyn_prop.set_integrator(pnt.RKF45)
        dyn_prop.set_integrator_params(
            pnt.IntegratorParams(max_iter=20, abstol=1e-10, reltol=1e-10)
        )
        dyn_prop.set_print_progress(False)
        self.dyn_gnss = dyn_prop

        # GNSS receiver parameters
        self.gnssr_param = GNSSReceiverParam()

        # create data directory if not exists
        gps_data_dir = self.basepath + "/gnss_meas"
        if not os.path.exists(gps_data_dir):
            os.makedirs(gps_data_dir)
            print(f"Created directory: {gps_data_dir}")

        self.gps_datadir = gps_data_dir

        print("Setup GNSS measurements directory complete.")

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
        P_tx_dict,
        freq="L1",
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
        P_tx = np.zeros((N_sc, N_gnss, N_t))

        # sun direction
        # rv_e2s_eci = pnt.get_body_pos_vel(t_tai, pnt.EARTH, pnt.SUN, pnt.ECI)
        t_tdb = pnt.convert_time(t_tai, pnt.TAI, pnt.TDB)
        rv_m2e_ci = pnt.get_body_pos_vel(t_tdb, pnt.MOON, pnt.EARTH, pnt.MOON_CI)

        self.rx_antenna = pnt.Antenna("moongpsr")

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

            print(f"    Computing visibility for spacecraft {i+1}/{N_sc}...")

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

                u_gps2sc = normalize(rv_m2sc_ci[i, :, :3] - rv_gnss_ci[j, :, :3])
                phi_gps2sc[i, j] = np.arccos(np.clip(np.sum(u_gps2sc * ez_gnss[j], axis=-1), -1, 1))
                theta_gps2sc[i, j] = np.arctan2(
                    np.sum(u_gps2sc * ey_gnss[j], axis=-1),
                    np.sum(u_gps2sc * ex_gnss[j], axis=-1),
                )

                u_sc2gps = -u_gps2sc
                u_sc2e = normalize(rv_m2e_ci[:, :3] - rv_m2sc_ci[i, :, :3])
                phi_sc2gps[i, j] = np.arccos(np.clip(np.sum(u_sc2e * u_sc2gps, axis=-1), -1, 1))

                G_tx[i, j] = gnss_antennas[prn_gnss[j]].compute_gain(
                    theta_gps2sc[i, j], phi_gps2sc[i, j]
                )
                G_rx[i, j] = self.rx_antenna.compute_gain(0, phi_sc2gps[i, j])

                G_tx[i, j][~vis_sc2gps[i, j]] = np.nan
                G_rx[i, j][~vis_sc2gps[i, j]] = np.nan

                P_tx[i, j] = P_tx_dict[prn_gnss[j]]

            vis_sc2sc[i, i, :] = False

        # Link budget -------------------------------------------------------
        CN0 = linkbudget_receiver(freq, dist_sc2gps, P_tx, G_tx, G_rx)

        return CN0, phi_gps2sc, theta_gps2sc, vis_sc2gps, G_tx, G_rx

    def compute_visibility(
        self, r1: np.ndarray, r2: np.ndarray, R_body: float, r_body: np.ndarray = None
    ) -> np.ndarray:

        Nt = r1.shape[0]  # number of time steps
        if r_body is None:
            r_body = np.zeros((Nt, 3))

        r1body_norm = np.linalg.norm(r1 - r_body, axis=1)
        r2body_norm = np.linalg.norm(r2 - r_body, axis=1)

        min_alt = R_body + 50e3  # minimum altitude for visibility in km
        visibility = np.ones(len(r1), dtype=bool)

        min_elev = np.deg2rad(-10)  # minimum elevation angle in radians

        r1_is_surface = np.any(r1body_norm < min_alt)
        r2_is_surface = np.any(r2body_norm < min_alt)

        # print(f"r1_is_surface: {r1_is_surface}, r2_is_surface: {r2_is_surface}")

        if r1_is_surface:
            # r1 is on surface -> judge by elevation angle
            r12 = r2 - r1  # N x 3
            r12_norm = np.linalg.norm(r12, axis=1)  # N
            r1_to_body = r_body - r1  # N x 3
            r1_to_body_norm = np.linalg.norm(r1_to_body, axis=1)  # N
            elev = np.zeros(len(r12_norm))  # elevation angle (N)

            for ti in range(Nt):
                tmp = (
                    np.dot(r12[ti, :], r1_to_body[ti, :]) / r12_norm[ti] / r1_to_body_norm[ti]
                )  # elevation angle (N)
                elev[ti] = np.arccos(np.clip(tmp, -1, 1)) - np.pi / 2

            # print(f"elev: {np.max(np.rad2deg(elev))}")
            visibility = elev > min_elev

        elif r2_is_surface:
            # r2 is on surface -> judge by elevation angle
            r21 = r1 - r2  # N x 3
            r21_norm = np.linalg.norm(r21, axis=1)  # N
            r2_to_body = r_body - r2
            r2_to_body_norm = np.linalg.norm(r2_to_body, axis=1)
            elev = np.zeros(Nt)
            for ti in range(len(r21_norm)):
                tmp = np.dot(r21[ti, :], r2_to_body[ti, :]) / r21_norm[ti] / r2_to_body_norm[ti]
                elev[ti] = np.arccos(np.clip(tmp, -1, 1)) - np.pi / 2

            visibility = elev > min_elev

        else:
            r = r2 - r1  # N x 3
            r_norm = np.linalg.norm(r, axis=1)
            r1body = r1 - r_body
            r1body_norm = np.linalg.norm(r1body, axis=1)
            dot = np.einsum("ij,ij->i", -r1body, r)
            theta1 = np.arccos(np.clip(dot / r_norm / r1body_norm, -1, 1))
            theta2 = np.arcsin(np.clip(R_body / r1body_norm, -1, 1))
            visibility[(theta1 < theta2) & (r_norm > r1body_norm)] = False

        return visibility

    def setup_gnss(
        self, gnss_consts=["GPS"], gps_datetime=None, overwrite=False, sp3_prop_method="interp"
    ):
        self.n_meas = 0
        t_tai = self.t_tai

        N_t = t_tai.shape[0]
        self.gnss_consts = gnss_consts
        self.gps_datetime = gps_datetime

        pnt.set_lupnt_epoch(0)

        # variables to store the results
        self.prn_gnss = {}
        self.rv_gnss_eci = {}
        self.rv_gnss_ci = {}
        self.rv_gnss_ecef = {}
        self.true_range = {}
        self.true_rangerate = {}
        self.true_carrier_phase = {}
        self.sigma_range = {}
        self.sigma_rangerate = {}
        self.sigma_carrier_phase = {}
        self.cn0 = {}
        self.vis_gnss = {}
        self.phi_gps2sc = {}
        self.theta_gps2sc = {}
        self.G_tx = {}
        self.G_rx = {}
        self.vis_sc2gps = {}
        self.N_gnss = {}
        self.antennas = {}
        self.P_tx = {}

        # Load SP3 files
        sp3l = SP3Loader(target_dt=gps_datetime, sim_t=t_tai[-1] - t_tai[0], dt_timesys=pnt.TAI)
        self.sp3l = sp3l
        self.brdc = BRDCLoader(
            target_dt=gps_datetime, sim_t=t_tai[-1] - t_tai[0], dt_timesys=pnt.TAI
        )
        sp3l_ref_epoch = datetime_to_tai(gps_datetime)

        sats = sp3l.sats
        print(f"Loaded {len(sats)} satellites from SP3 file.")

        for gnss_const in gnss_consts:
            print(" ")
            print(f"Setting up {gnss_const} ===========================")

            # find satellites that starts with 'G' for GPS, 'E' for Galileo
            if gnss_const == "GPS":
                sats_gnss = [s for s in sats if s.startswith("G")]
                prns = [int(s[1:]) for s in sats_gnss]
                consts = [s[0] for s in sats_gnss]
            elif gnss_const == "GALILEO":
                sats_gnss = [s for s in sats if s.startswith("E")]
                prns = [int(s[1:]) for s in sats_gnss]
                consts = [s[0] for s in sats_gnss]
            elif gnss_const == "QZSS":
                sats_gnss = [s for s in sats if s.startswith("J")]
                prns = [int(s[1:]) for s in sats_gnss]
                consts = [s[0] for s in sats_gnss]
            else:
                raise ValueError(f"Unsupported GNSS constellation: {gnss_const}")

            N_gnss = len(sats_gnss)
            rv_gnss_eci = np.zeros((N_gnss, N_t, 6))
            rv_gnss_ci = np.zeros((N_gnss, N_t, 6))
            rv_gnss_ecef = np.zeros((N_gnss, N_t, 6))

            print("propagating {} orbits...".format(N_gnss))
            print(" prns:", prns)

            sp3l_ref_epoch = np.median(sp3l.epochs)

            gps_orbit_file = (
                self.gps_datadir
                + "/{}_orbits_date_{}_{}_{}_{}_Nt_{}_dt_{}.npz".format(
                    gnss_const.lower(),
                    gps_datetime.year,
                    gps_datetime.month,
                    gps_datetime.day,
                    gps_datetime.hour,
                    N_t,
                    int(t_tai[1] - t_tai[0]),
                )
            )

            if os.path.exists(gps_orbit_file) and not overwrite:
                print(f"Loading {gnss_const} orbits from {gps_orbit_file}...")
                data = np.load(gps_orbit_file)
                rv_gnss_eci = data["rv_gnss_eci"]
                rv_gnss_ci = data["rv_gnss_ci"]
                rv_gnss_ecef = data["rv_gnss_ecef"]
                prns = data["prns"]
            else:
                for i in range(N_gnss):
                    print(f"  satellite {i+1}/{N_gnss}")
                    # use the sp3 epochs for earth orbits
                    if sp3_prop_method == "propagate":
                        # proagate from the reference epoch using two-body dynamics
                        rv0_gps_eci = sp3l.get_posvel(
                            consts[i], prns[i], sp3l_ref_epoch, out_frame=pnt.ECI
                        )
                        t_prop = t_tai - sp3l_ref_epoch
                        pnt.set_lupnt_epoch(sp3l_ref_epoch)
                        rv_gnss_eci[i] = self.dyn_gnss.propagate(rv0_gps_eci, t_prop)
                        pnt.set_lupnt_epoch(0)
                        t_tdb = pnt.convert_time(t_tai, pnt.TAI, pnt.TDB)
                        rv_gnss_ecef[i] = pnt.convert_frame(
                            t_tdb, rv_gnss_eci[i], pnt.ECI, pnt.ECEF
                        )
                        rv_gnss_ci[i] = pnt.convert_frame(
                            t_tdb, rv_gnss_eci[i], pnt.ECI, pnt.MOON_CI
                        )
                    elif sp3_prop_method == "interp":
                        # a more accurate retrieval from SP3 file at each time step
                        rv_gnss_ecef[i] = sp3l.get_posvel(
                            consts[i], prns[i], t_tai, out_frame=pnt.ECEF
                        )
                        t_tdb = pnt.convert_time(t_tai, pnt.TAI, pnt.TDB)
                        rv_gnss_ci[i] = pnt.convert_frame(
                            t_tdb, rv_gnss_ecef[i], pnt.ECEF, pnt.MOON_CI
                        )
                        rv_gnss_eci[i] = pnt.convert_frame(
                            t_tdb, rv_gnss_ecef[i], pnt.ECEF, pnt.ECI
                        )
                    else:
                        raise ValueError("sp3_prop_method must be 'propagate' or 'interp'")

                # save the orbits to a file
                print(f"Saving {gnss_const} orbits to {gps_orbit_file}...")
                np.savez(
                    gps_orbit_file,
                    rv_gnss_eci=rv_gnss_eci,
                    rv_gnss_ci=rv_gnss_ci,
                    rv_gnss_ecef=rv_gnss_ecef,
                    prns=prns,
                )

            # antennas
            print(f"Setting up {gnss_const} antennas...")
            antennas = {}
            P_tx = {}

            if gnss_const == "GPS":
                df = pd.read_csv(pnt.get_file_path("gps_table.csv"))
                print("df_size:", df.shape)
                gps_antenna_names = df.set_index("PRN")["LM_File"].to_dict()
                block_names = df.set_index("PRN")["blockName"].to_dict()
                # for gps that LM file is not available, use the ACE_file instead (IIF only has ACE file)
                for prn in gps_antenna_names.keys():
                    if type(gps_antenna_names[prn]) == float:
                        gps_antenna_names[prn] = df.loc[prn - 1, "ACE_File"]

                print(gps_antenna_names)
                print(block_names)

                freqs = ["L1", "L5"]
                for freq in freqs:
                    antennas[freq] = {k: pnt.Antenna(v) for k, v in gps_antenna_names.items()}
                    P_tx[freq] = {k: block_power(block_names[k], freq) for k in prns}

            elif gnss_const == "GALILEO":
                freqs = ["E1", "E5a"]
                for freq in freqs:
                    antennas[freq] = {k: pnt.Antenna("Galileo_{}".format(freq)) for k in prns}
                    P_tx[freq] = {k: 14.0 for k in prns}  # Fixed  value for Galileo

            elif gnss_const == "QZSS":
                qzss_names = ["1R", "02", "03", "04", "05", "06", "07"]
                freqs = ["L1", "L5"]
                for freq in freqs:
                    if freq == "L1":
                        offset = 0.0
                    elif freq == "L5":
                        offset = 3.0
                    antennas[freq] = {
                        k: pnt.Antenna("QZSS_" + qzss_names[k - 1] + "_" + freq) for k in prns
                    }

                    # https://navi.ion.org/content/71/2/navi.641
                    P_tx[freq] = {k: 14.1 + offset for k in prns}  # Fixed value for QZSS

            else:
                raise ValueError(f"Unsupported GNSS constellation: {gnss_const}")

            # store the results
            print(f"Storing {gnss_const} data...")
            self.prn_gnss[gnss_const] = prns
            self.rv_gnss_eci[gnss_const] = rv_gnss_eci
            self.rv_gnss_ci[gnss_const] = rv_gnss_ci
            self.rv_gnss_ecef[gnss_const] = rv_gnss_ecef
            self.N_gnss[gnss_const] = N_gnss
            self.antennas[gnss_const] = antennas
            self.P_tx[gnss_const] = P_tx

            print(" ")

    def setup_measurements(self, cn0_threshold=15.0, overwrite=False):
        # sun direction
        t_tdb = pnt.convert_time(self.t_tai, pnt.TAI, pnt.TDB)
        rv_e2s_eci = pnt.get_body_pos_vel(t_tdb, pnt.EARTH, pnt.SUN, pnt.ECI)

        # constants
        N_t = self.t_tai.shape[0]
        dt = self.t_tai[1] - self.t_tai[0]
        gps_datetime = self.gps_datetime
        N_sc = self.rv_m2sc_ci.shape[0]  # number of spacecraft
        print(f"Setting up measurements for {N_sc} spacecraft and {N_t} time steps =========")

        for gnss_const in self.gnss_consts:
            if gnss_const == "GPS":
                # freqs = ["L1", "L2", "L5"]
                freqs = ["L1", "L5"]
                # https://www.navcen.uscg.gov/gps-constellation
                # P_tx = 14.3
            elif gnss_const == "GALILEO":
                # freqs = ["E1", "E5a", "E5b", "E6"]
                freqs = ["E1", "E5a"]
                # https://www.gsc-europa.eu/system-service-status/constellation-information
                # P_tx = 14.1
            elif gnss_const == "QZSS":
                freqs = ["L1", "L5"]
                # P_tx = 14.1
            else:
                raise ValueError(f"Unsupported GNSS constellation: {gnss_const}")

            fault_prns = get_fault_prns(gnss_const, use_faults=self.consider_faults)

            fault_idxs = []
            for prn in fault_prns:
                # print(f"Checking for faulty PRN {prn} in {gnss_const}...")
                if prn in self.prn_gnss[gnss_const]:
                    prns = np.atleast_1d(self.prn_gnss[gnss_const])
                    match = np.where(prns == prn)[0]
                    if match.size == 0:
                        raise ValueError(
                            f"PRN {prn} not found in prn_gnss[{gnss_const}] (shape={np.shape(self.prn_gnss[gnss_const])})"
                        )
                    fault_idxs.append(int(match[0]))

            print(f"Setting up {gnss_const} measurements...")
            print(f"Faulty PRNs for {gnss_const}: {fault_prns} at indices {fault_idxs}")
            fault_idxs = []

            N_gnss = self.N_gnss[gnss_const]
            rv_gnss_ci = self.rv_gnss_ci[gnss_const]
            rv_gnss_eci = self.rv_gnss_eci[gnss_const]
            prns = self.prn_gnss[gnss_const]
            antennas = self.antennas[gnss_const]
            P_tx = self.P_tx[gnss_const]

            # attitde computation
            e_gnss2e = normalize(-rv_gnss_eci[:, :, :3])
            e_gnss2s = normalize(rv_e2s_eci[None, :, :3] - rv_gnss_eci[:, :, :3])
            ez_gnss = e_gnss2e
            ey_gnss = cross_norm(e_gnss2e, e_gnss2s)
            ex_gnss = cross_norm(ey_gnss, ez_gnss)

            # Compute the CN0
            self.cn0[gnss_const] = {}
            self.vis_gnss[gnss_const] = {}
            self.phi_gps2sc[gnss_const] = {}
            self.theta_gps2sc[gnss_const] = {}
            self.G_tx[gnss_const] = {}
            self.G_rx[gnss_const] = {}
            self.vis_sc2gps[gnss_const] = {}
            self.sigma_range[gnss_const] = {}
            self.sigma_rangerate[gnss_const] = {}
            self.sigma_carrier_phase[gnss_const] = {}

            for freq in freqs:

                # store the results to a file
                gps_orbit_file = (
                    self.gps_datadir
                    + "/{}_{}_cn0_{}_meas_date_{}_{}_{}_{}_Nt_{}_dt_{}.npz".format(
                        gnss_const.lower(),
                        freq,
                        cn0_threshold,
                        gps_datetime.year,
                        gps_datetime.month,
                        gps_datetime.day,
                        gps_datetime.hour,
                        N_t,
                        int(dt),
                    )
                )

                if os.path.exists(gps_orbit_file) and not overwrite:
                    file_exists = True
                else:
                    file_exists = False

                print(f"  Frequency {freq}... {file_exists}")

                if file_exists:
                    print(f"Loading {gnss_const} {freq} measurements from {gps_orbit_file}...")
                    data = np.load(gps_orbit_file)
                    self.cn0[gnss_const][freq] = data["cn0"]
                    self.vis_gnss[gnss_const][freq] = data["vis_gnss"]
                    self.phi_gps2sc[gnss_const][freq] = data["phi_gps2sc"]
                    self.theta_gps2sc[gnss_const][freq] = data["theta_gps2sc"]
                    self.G_tx[gnss_const][freq] = data["G_tx"]
                    self.G_rx[gnss_const][freq] = data["G_rx"]
                    self.vis_sc2gps[gnss_const][freq] = data["vis_sc2gps"]
                    self.sigma_range[gnss_const][freq] = data["sigma_range"]
                    self.sigma_rangerate[gnss_const][freq] = data["sigma_rangerate"]
                    self.sigma_carrier_phase[gnss_const][freq] = data["sigma_carrier_phase"]

                else:
                    CN0, phi_gps2sc, theta_gps2sc, vis_sc2gps, G_tx, G_rx = self.compute_cn0(
                        self.t_tai,
                        self.rv_m2sc_ci,
                        N_gnss,
                        rv_gnss_ci,
                        ex_gnss,
                        ey_gnss,
                        ez_gnss,
                        antennas[freq],
                        prns,
                        P_tx[freq],
                        freq,
                    )

                    vis_gnss = (
                        CN0 > cn0_threshold
                    )  # N_sc x N_gnss x N_t visibility mask based on CN0 threshold

                    # set unhealthy satellites to NaN
                    vis_gnss[:, fault_idxs, :] = False

                    # compute the true range, range rate, and carrier phase
                    sigma_range = np.zeros((N_sc, N_gnss, N_t))
                    sigma_rangerate = np.zeros((N_sc, N_gnss, N_t))
                    sigma_carrier_phase = np.zeros((N_sc, N_gnss, N_t))
                    sigma_range[vis_gnss] = self.compute_gnss_pseudorange_noise(
                        CN0[vis_gnss], self.gnssr_param, freq
                    )
                    sigma_rangerate[vis_gnss] = self.compute_gnss_pseudorangerate_noise(
                        CN0[vis_gnss], freq
                    )
                    sigma_carrier_phase[vis_gnss] = self.compute_gnss_carrier_phase_noise(
                        CN0[vis_gnss], freq
                    )

                    # storage
                    self.cn0[gnss_const][freq] = CN0
                    self.vis_gnss[gnss_const][freq] = vis_gnss

                    self.phi_gps2sc[gnss_const][freq] = phi_gps2sc
                    self.theta_gps2sc[gnss_const][freq] = theta_gps2sc
                    self.G_tx[gnss_const][freq] = G_tx
                    self.G_rx[gnss_const][freq] = G_rx
                    self.vis_sc2gps[gnss_const][freq] = vis_sc2gps

                    self.sigma_range[gnss_const][freq] = sigma_range
                    self.sigma_rangerate[gnss_const][freq] = sigma_rangerate
                    self.sigma_carrier_phase[gnss_const][freq] = sigma_carrier_phase

                    # save the results to a file
                    print(f"Saving {gnss_const} {freq} measurements to {gps_orbit_file}...")
                    np.savez(
                        gps_orbit_file,
                        cn0=CN0,
                        vis_gnss=vis_gnss,
                        phi_gps2sc=phi_gps2sc,
                        theta_gps2sc=theta_gps2sc,
                        G_tx=G_tx,
                        G_rx=G_rx,
                        vis_sc2gps=vis_sc2gps,
                        sigma_range=sigma_range,
                        sigma_rangerate=sigma_rangerate,
                        sigma_carrier_phase=sigma_carrier_phase,
                    )

            # true range, range rate, and carrier phase
            # rv_m2sc_ci = [N_sc, N_t, 6] spacecraft position in Moon CI frame
            # rv_gnss_ci = [N_gnss, N_t, 6] GNSS satellite position in Moon CI frame
            print("  True range, range rate, and carrier phase...")
            true_range = np.zeros((N_sc, N_gnss, N_t))
            true_rangerate = np.zeros((N_sc, N_gnss, N_t))
            true_carrier_phase = np.zeros((N_sc, N_gnss, N_t))
            for i in range(N_sc):
                for j in range(N_gnss):
                    # Todo: solve light time correction
                    # true range
                    true_range[i, j] = np.linalg.norm(
                        self.rv_m2sc_ci[i, :, :3] - rv_gnss_ci[j, :, :3], axis=1
                    )

                    # true range rate
                    true_rangerate[i, j] = (
                        np.sum(
                            (self.rv_m2sc_ci[i, :, 3:] - rv_gnss_ci[j, :, 3:])
                            * (self.rv_m2sc_ci[i, :, :3] - rv_gnss_ci[j, :, :3]),
                            axis=1,
                        )
                        / true_range[i, j]
                    )

                    # true carrier phase
                    true_carrier_phase[i, j] = true_range[i, j]

            print(" ")

    def compute_gnss_pseudorange_noise(self, CN0_dB, gnssr_param, freq):
        CN0 = 10 ** (CN0_dB / 10)
        Bn = gnssr_param.Bn
        Rc = gnss_rc_map[freq]  # Chip rate in Hz
        Bfe = gnssr_param.b * Rc
        T = gnssr_param.T
        D = gnssr_param.D
        Tc = 1.0 / Rc
        C = pnt.C  # Speed of light in m/s

        sigma = np.zeros_like(CN0)

        case1 = D >= (np.pi * Rc / Bfe)
        case2 = (D > (Rc / Bfe)) & (~case1)
        case3 = ~case1 & ~case2

        # Case 1: Wide spacing discriminator
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
        # https://theses.eurasip.org/media/theses/documents/padma-bolla-advanced-tracking-loop-architectures-for-multi-frequency-gnss-receiver.pdf?utm_source=chatgpt.com
        gnssr_param = self.gnssr_param
        f = gnss_freq_map[freq]
        lambda_ = pnt.C / f  # Wavelength of the GNSS signal

        CN0 = 10 ** (CN0_dB / 10)
        F = 2
        Bf = gnssr_param.Bf
        T = gnssr_param.T

        return lambda_ / (2 * np.pi * T) * np.sqrt(4 * Bf / CN0 * (1 + 1.0 / (T * CN0)))

    def compute_gnss_carrier_phase_noise(self, CN0_dB, freq):
        CN0 = 10 ** (CN0_dB / 10)
        f = gnss_freq_map[freq]
        lambda_ = pnt.C / f  # Wavelength of the GNSS signal

        gnssr_param = self.gnssr_param
        Bp = gnssr_param.Bp
        T = gnssr_param.T

        return lambda_ / (2 * np.pi) * np.sqrt(Bp / CN0 * (1.0 + 1.0 / (2 * T * CN0)))

    def lt_correction(self, gnss_const, sat_idx, gps_idx, tidx):
        """
        Compute the light time correction for the given time index and receiver position.
        """
        # compute the distance from the receiver to the GNSS satellite
        rv_tx_eci0 = self.rv_gnss_eci[gnss_const][gps_idx, tidx, :6]  # [6]
        pos_tx_eci = rv_tx_eci0[:3]  # [3]

        pos_rx_ecef = self.rv_e2sc_ecef[sat_idx, tidx, :3]  # [3]
        t_rx_tdb = pnt.convert_time(self.t_tai[tidx], pnt.TAI, pnt.TDB)
        pos_rx_eci = pnt.convert_frame(t_rx_tdb, pos_rx_ecef, pnt.ECEF, pnt.ECI, rotate_only=False)

        # distance
        dist = np.linalg.norm(pos_rx_eci - pos_tx_eci)  # [N_sc, N_gnss]

        # compute the light time correction
        n_iter = 3

        for i in range(n_iter):
            lt = dist / pnt.C
            rv_tx_eci = self.dyn_gnss.propagate(
                rv_tx_eci0, self.t_tai[tidx], self.t_tai[tidx] - lt, np.zeros(6)
            )
            pos_tx_eci = rv_tx_eci[:3]  # [3]
            new_dist = np.linalg.norm(pos_tx_eci - pos_rx_eci)  # [N_sc, N_gnss]

            if np.linalg.norm(new_dist - dist) < 1e-3:  # under 1 mm
                break

            dist = new_dist

        # convert to ECEF frame
        t_tx_tdb = pnt.convert_time(self.t_tai[tidx] - lt, pnt.TAI, pnt.TDB)
        rv_tx_ecef = pnt.convert_frame(t_tx_tdb, rv_tx_eci, pnt.ECI, pnt.ECEF, rotate_only=False)

        return rv_tx_ecef, lt

    def plot_gnss_orbit(self, savefig=False, filename=None, plot_inv=10, plot_moon_dir=True):
        tickvals = np.arange(-30, 31, 10) * 1e3
        fig = go.Figure()

        pnt.plot.plot_body(fig, pnt.EARTH, size_factor=5)

        orbit_colors = {
            "GPS": "blue",
            "GALILEO": "orange",
            "QZSS": "green",
        }

        # gnss
        for gnss_const in self.gnss_consts:
            print("Plotting {} orbits...".format(gnss_const))
            rv_eci = self.rv_gnss_eci[gnss_const]
            tidx = np.arange(0, rv_eci.shape[1], plot_inv)

            pnt.plot.plot_orbits(fig, rv_eci[:, tidx, :], t=0, color=orbit_colors[gnss_const])

        # eartg
        if plot_moon_dir:
            print("Plotting Moon direction...")
            # moon direction
            t0_tdb = pnt.convert_time(self.t_tai[0], pnt.TAI, pnt.TDB)
            rv_e2m_eci = pnt.get_body_pos_vel(t0_tdb, pnt.EARTH, pnt.MOON, pnt.ECI)

            rv_e2m_eci_line = np.zeros((100, 3))
            scale_factor = 10  # shorten the line segment
            for i in range(100):
                rv_e2m_eci_line[i, :] = i * rv_e2m_eci[:3] / 100 / scale_factor

            # a line segemnt from Earth to Moon
            pnt.plot.plot_orbits(fig, rv_e2m_eci_line[None, :, :3], t=0, color="gray")

        pnt.plot.set_view(fig, -80, 20, 2.5)

        fig.update_layout(showlegend=True, width=400, height=400)

        if savefig:
            print("Saving figure...")
            if filename is None:
                filename = self.basepath + "/figures/gnss_orbit.pdf"
            fig.write_image(filename)
            print(f"Figure saved to {filename}")

        return fig

    def plot_num_tracked_sats(self, sat_labels=None, savefig=False, filename=None):

        n_const = len(self.gnss_consts)
        n_sats = self.rv_m2sc_ci.shape[0]
        tspan = self.t_tai - self.t_tai[0]

        fig, axes = plt.subplots(n_const, n_sats, figsize=(16, 8))

        for i, gnss_const in enumerate(self.gnss_consts):

            if gnss_const == "GPS":
                signal = "L1"
            elif gnss_const == "GALILEO":
                signal = "E1"

            for j in range(n_sats):
                ax = axes[i, j]
                ax.plot(
                    tspan / 3600,
                    np.sum(self.vis_gnss[gnss_const][signal][j], axis=0),
                    label=gnss_const,
                )
                ax.set_title(f"{gnss_const} - {sat_labels[j] if sat_labels else 'Sat ' + str(j+1)}")
                ax.set_xlabel("Time (TAI)")
                ax.set_ylabel("Visibility")
                ax.grid()
                ax.legend()

        plt.tight_layout()

        if savefig:
            if filename is None:
                filename = self.basepath + "/figures/num_tracked_sats.pdf"
            plt.savefig(filename)
            print(f"Figure saved to {filename}")

    def plot_histogram(self, cn0_threshold=20.0):
        """
        Plot the histogram of the number of tracked satellites.
        """
        metrics = [
            self.G_tx,  # G_tx
            self.cn0,  # CN0
            self.sigma_range,  # sigma_range
            self.phi_gps2sc,  # phi_gps2sc
        ]

        bins = [
            np.linspace(-30, 20, 50),  # G_tx
            np.linspace(0, 40, 50),  # CN0
            np.linspace(0, 20, 50),  # sigma_range
            np.linspace(0, 90, 50),  # phi_gps2sc
        ]
        xlabels = [
            "G_tx (dB)",  # G_tx
            "CN0 (dB-Hz)",  # CN0
            "Sigma Range (m)",  # sigma_range
            "Boresite Angle (deg)",  # phi_gps2sc
        ]

        scales = [1, 1, 1, 180 / np.pi]

        for i, metric in enumerate(metrics):
            fig, ax = plt.subplots(1, 1, figsize=(6, 4))
            colors = [
                "blue",
                "orange",
                "green",
                "red",
                "purple",
                "brown",
                "pink",
                "gray",
                "cyan",
                "magenta",
            ]
            labels = []
            plot_y = []

            for j, gnss_const in enumerate(self.gnss_consts):
                print(f"Plotting {gnss_const} {xlabels[i]}...")
                if gnss_const == "GPS":
                    signals = ["L1"]
                elif gnss_const == "GALILEO":
                    signals = ["E1"]

                for signal in signals:
                    vis = self.cn0[gnss_const][signal] > cn0_threshold
                    tmp = (
                        metric[gnss_const][signal][vis > 0].flatten() * scales[i]
                    )  # scale the data
                    plot_y.append(tmp)
                    labels.append(f"{gnss_const} {signal}")

            if i == 3:  # boresite angle
                for j, yconst in enumerate(plot_y):
                    gnss_const = self.gnss_consts[j]
                    min_y = np.min(yconst)
                    max_y = np.max(yconst)
                    print(f"  {gnss_const} boresite angle: min={min_y:.2f}, max={max_y:.2f}")

            ax.hist(
                plot_y,
                bins=bins[i],
                density=True,
                histtype="bar",
                label=labels,
                color=colors[: len(labels)],
                alpha=0.7,
            )
            ax.set_xlabel(xlabels[i])
            ax.set_ylabel("Density")
            ax.grid()
            ax.legend()

        plt.tight_layout()

    def plot_scatter(self, cn0_threshold=20.0):
        metrics = [
            self.G_tx,  # G_tx
            self.cn0,  # CN0
            self.sigma_range,  # sigma_range
            self.phi_gps2sc,  # phi_gps2sc
            self.vis_sc2gps,  # vis_sc2gps
        ]

        metric_labels = [
            "G_tx (dB)",  # G_tx
            "CN0 (dB-Hz)",  # CN0
            "Sigma Range (m)",  # sigma_range
            "Boresite Angle (deg)",  # phi_gps2sc
            "Visibility (True/False)",  # vis_sc2gps
        ]

        scales = [1, 1, 1000, 180 / np.pi, 1]

        plot_y = {}
        labels = {}

        for i, metric in enumerate(metrics):
            colors = [
                "blue",
                "orange",
                "green",
                "red",
                "purple",
                "brown",
                "pink",
                "gray",
                "cyan",
                "magenta",
            ]
            labels[i] = []
            plot_y[i] = []

            for j, gnss_const in enumerate(self.gnss_consts):
                if gnss_const == "GPS":
                    signals = ["L1"]
                elif gnss_const == "GALILEO":
                    signals = ["E1"]

                for signal in signals:
                    vis = self.cn0[gnss_const][signal] > cn0_threshold

                    tmp = (
                        metric[gnss_const][signal][vis > 0].flatten() * scales[i]
                    )  # scale the data
                    plot_y[i].append(tmp[::3])
                    labels[i].append(f"{gnss_const} {signal}")

                    max_y = np.max(tmp)
                    min_y = np.min(tmp)
                    print(f"{gnss_const} {metric_labels[i]}: min={min_y:.2f}, max={max_y:.2f}")

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        axes = axes.flatten()

        # boresite angle vs cn0 ----------------------------------
        ax = axes[0]
        for j in range(len(plot_y[3])):  # for each constellation
            ax.scatter(
                plot_y[3][j],  # boresite angle
                plot_y[1][j],  # CN0
                c=colors[j],
                alpha=0.5,
                label=labels[1][j],
            )
        ax.set_xlabel("Boresite Angle (deg)")
        ax.set_ylabel("CN0 (dB-Hz)")
        ax.grid()
        ax.legend()

        # Boresite vs G_tx ----------------------------------
        ax = axes[1]
        for j in range(len(plot_y[0])):  # for each constellation
            ax.scatter(
                plot_y[3][j],  # boresite angle
                plot_y[0][j],  # G_tx
                c=colors[j],
                alpha=0.5,
                label=labels[1][j],
            )
        ax.set_xlabel("Boresite Angle (deg)")
        ax.set_ylabel("G_tx (dB)")
        ax.grid()
        ax.legend()
