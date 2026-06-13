import pylupnt as pnt
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import os
from tqdm import tqdm


class LunanetSatAntenna:
    def __init__(self, P_tx, coe):
        self.setup_antenna(P_tx, coe)

    def setup_antenna(self, P_tx, coe):
        basepath = pnt.get_basepath()
        csv_file = os.path.join(basepath, "antenna", "LGPS", "lunanet_esa.csv")
        df = pd.read_csv(csv_file, header=None)
        self.boresite_angle = df.iloc[:, 0].values  # elevation angles in degrees
        self.boresite_angle[0] = 0.0  # set the first element to 0 degrees
        self.gains = df.iloc[:, 1].values  # antenna gains in dB

        # compute scaling factor
        R_MOON = pnt.R_MOON
        a = coe[0]
        e = coe[1]
        # print("antenna: a = {}, e = {}".format(a, e))
        a_nom = 9750e3  # from ESA design
        e_nom = 0.6383  # from ESA design
        nominal_angle = np.arctan(R_MOON / a_nom / (1 + e_nom))
        angle = np.arctan(R_MOON / a / (1 + e))  # angle at the periapsis

        self.scale_phi = angle / nominal_angle

        self.interp = CubicSpline(self.boresite_angle, self.gains)
        self.P_tx = P_tx  # default value

    def compute_eirp(self, phi_rad):
        phi_deg = np.rad2deg(phi_rad)
        idx = np.where((phi_deg >= self.boresite_angle[0]) & (phi_deg <= self.boresite_angle[-1]))[
            0
        ]
        phi_deg[idx] = phi_deg[idx] / self.scale_phi  # scale to the nominal design
        gain = self.interp(phi_deg[idx])

        eirp = -np.inf * np.ones_like(phi_rad)
        eirp[idx] = gain + self.P_tx  # EIRP in dBW

        return eirp

    def plot_pattern(self):
        """
        Plot the antenna gain pattern.
        """
        fig = plt.figure(figsize=(6, 4))
        eps = 1e-6  # small value to avoid division by zero
        phi_rads = np.linspace(0 + eps, np.pi / 2 - eps, 500)  # from 0 to 90 degrees
        gains = self.compute_eirp(phi_rads)
        plt.plot(np.rad2deg(phi_rads), gains, "bo-", label="Antenna Gain Pattern")
        plt.xlabel("Elevation Angle (degrees)")
        plt.ylabel("EIRP (dBW)")
        plt.title("Lunanet Antenna Gain Pattern")
        plt.grid()
        plt.xlim(0, 90)
        plt.ylim(-50, 50)


class LunanetReceiverParam:
    def __init__(self):
        # GNSS Receiver chip parameter
        # https://www.mdpi.com/1424-8220/16/3/347?utm_source=chatgpt.com
        self.Bp = 5.0  # Carrier loop noise bandwidth [Hz]
        self.T = 20e-3  # Tracking loop integration time [s]
        self.b = 2.0  # normalized bandwidth [Hz]
        self.Bn = 1.0  # Code loop noise bandwidth [Hz]
        self.Bf = 10.0  # Frequency loop noise bandwidth [Hz]
        self.D = 1.0  # Early-to-late correlator spacing (chips)

        self.lambda_c = 58.61e-3  # [km] code wavelength
        C = 299792.458  # [km/s] speed of light
        self.Rc = C / self.lambda_c  # chip rate   C/Rc = lambda_c
        self.freq_Hz = 2491.005e6  # [Hz] frequency


def gps_patch_gain_dbi(
    elev_rad,
    G_zenith=3.0,  # dBi at zenith (typical +3 to +6)
    G_horizon=-5.0,  # dBi at horizon (typical 0 to -5)
    q=2.0,  # roll-off exponent (2=cos^2; higher = sharper)
    back_lobe=-15.0,  # dBi below horizon
):
    """
    Elevation-dependent receive gain for a hemispherical RHCP GPS patch.
    elev_rad : elevation angle(s) in radians (0 at horizon, +pi/2 at zenith)
    Returns G_rx in dBi (numpy array or float).
    """
    elev = np.asarray(elev_rad, dtype=float)
    # Off-boresight from zenith
    phi = (np.pi / 2.0) - elev

    # Cosine roll-off toward the horizon, clamped to 0 (no negative cos)
    c = np.cos(np.clip(phi, 0.0, np.pi / 2.0))
    shape = c**q  # 1 at zenith, 0 at horizon

    # Interpolate between zenith and horizon gains
    G = G_horizon + (G_zenith - G_horizon) * shape

    # Apply back-lobe value for below-horizon rays
    G = np.where(elev >= 0.0, G, back_lobe)
    return G


def compute_lunanet_cn0(
    t_tai,
    x_orb,
    r_user2sat,
    elev_user,
    user2sat_norm,
    lunanet_antennas,
    min_elev_deg=10.0,
    cn0_thresh=30.0,
):
    """
    Compute CN0 values and determine coverage from a satellite to surface user.

    Parameters:
        t

    Returns:
        cn0 : np.ndarray [n_user, n_sat, lent]
        coverage : np.ndarray [n_user, n_sat, lent]
    """

    # Geometry --------------------------------------------------------
    n_user, n_sat, lent = r_user2sat.shape[:3]

    # Satellite boresight assumed toward Moon center
    e_sat2moon = -x_orb[:, :, :3]
    e_sat2moon = e_sat2moon / np.linalg.norm(
        e_sat2moon, axis=-1, keepdims=True
    )  # normalize (n_sat, lent, 3)

    # Broadcast satellite boresight to user shape (n_user, n_sat, lent, 3)
    ez_sat = np.tile(e_sat2moon[np.newaxis, :, :, :], (n_user, 1, 1, 1))

    # Unit vector from satellite to user
    u_sat2user = -r_user2sat / user2sat_norm[..., np.newaxis]  # (n_user, n_sat, lent, 3)

    # Angle from satellite boresight to user direction
    phi_sat = np.arccos(np.clip(np.sum(u_sat2user * ez_sat, axis=-1), -1.0, 1.0))  # [rad]

    # Receiver model --------------------------------------------------
    T_ant = 100.0  # [K] Antenna temperature
    cable_loss_before_LNA = 1.0  # [dB]
    LNA_NF = 3.0  # [dB]
    LNA_gain = 30.0  # [dB]
    cable_loss_after_LNA = 10.0  # [dB]
    T0 = 290.0  # [K]

    # Convert to linear units
    LNA_F = 10 ** (LNA_NF / 10)
    L1 = 10 ** (cable_loss_before_LNA / 10)
    L2 = 10 ** (cable_loss_after_LNA / 10)
    G_LNA = 10 ** (LNA_gain / 10)

    # Effective noise temperature
    T_LNA = T0 * (LNA_F - 1)
    T_cable1 = T0 * (L1 - 1)
    T_anteff = L1 * T_ant
    T_cable2 = T0 * (L2 - 1)
    Teff = T_anteff + T_cable1 + T_LNA / L1 + T_cable2 / (L1 * G_LNA)  # [K]

    # Link budget -----------------------------------------------------
    freq_Hz = 2491.005e6
    c = getattr(pnt, "C", 299792458.0)  # fallback if pnt.C is not defined

    # Free-space loss [dB]
    L_fs = 20 * np.log10((4 * np.pi * user2sat_norm * freq_Hz) / c)

    # Boltzmann constant in dB
    kb_db = 10 * np.log10(1.38e-23)

    # EIRP computation from antenna model
    EIRP_tx = np.zeros((n_user, n_sat, lent))
    cable_loss_tx = 1.0  # [dB]

    for i in range(n_sat):
        phi_i = phi_sat[:, i, :].flatten()
        eirp_i = lunanet_antennas[i].compute_eirp(phi_i).reshape(n_user, lent)
        EIRP_tx[:, i, :] = eirp_i

    # Parameters
    G0_rx = 6.0  # [dB] maximum gain (on-axis)
    rx_hpbw = np.deg2rad(35.0)  # [rad] half-power beamwidth
    phi_user = np.pi / 2 - elev_user  # [rad] off-boresight angle (0 at zenith)

    type = "parabolic"
    # type = "patch"

    if type == "patch":
        G_rx = gps_patch_gain_dbi(elev_user)
    else:
        G_rx = G0_rx - 3.0 * (phi_user / rx_hpbw) ** 2

    # Optional: apply a sidelobe floor (e.g., -10 dB below peak)
    G_rx = np.maximum(G_rx, G0_rx - 20.0)

    # Compute CN0 [dB-Hz]
    cn0 = EIRP_tx - cable_loss_tx - L_fs - kb_db + G_rx - 10 * np.log10(Teff)

    # Apply coverage masks
    coverage_elev = elev_user > np.deg2rad(min_elev_deg)
    coverage_cn0 = cn0 > cn0_thresh
    coverage = coverage_elev & coverage_cn0

    return cn0, coverage
