"""Lunar GNSS simulation configuration helper (``pylupnt.config``).

Build a ready-to-run :class:`GNSSMeasurements` pipeline for a lunar GNSS
receiver without wiring every option by hand.  A single flag selects the transmit
model:

``mode="native"``
    LuPNT's default per-SVN **3-D** antenna pattern for each satellite (Lockheed-Martin
    panel patterns for GPS IIR/IIR-M/III, GPS ACE for IIA/IIF, GRAP for Galileo),
    evaluated at the modeled yaw-steered azimuth (``FULL_2D``), with LuPNT's default
    (higher) transmit power from the LuGRE payload analysis.

``mode="reference"``
    A reference-simulator setting: the azimuth-averaged **ACE block-average**
    patterns (a "2-D" off-boresight pattern) with the lower per-block transmit power that
    reproduces the reference C/N0.

The GNSS constellation is seeded from the latest precise SP3 ephemeris and propagated with
full dynamics (Earth 8x8 + Sun + Moon, RK8) to the requested epochs -- so **future** epochs
(with no SP3/broadcast product) can be simulated.  The receiver uses the stored ``moongpsr``
antenna and the reference link-budget settings.

Example
-------
>>> import numpy as np, pylupnt as pnt
>>> cfg = pnt.config.LunarGnssConfig(mode="reference")
>>> epochs_tai = t0_tai + np.arange(0.0, 360.0, 6.0)          # 6 min @ 6 s
>>> gm = cfg.build_measurements(epochs_tai)
>>> channels = gm.build_channels(float(epochs_tai[0]), rx_state_eci)   # per epoch
"""

from __future__ import annotations

import csv

import numpy as np

import pylupnt as pnt

__all__ = ["LunarGnssConfig", "lunar_gnss_measurements", "LUNAR_RX_PARAMS"]

# --- Lunar GNSS receiver link-budget settings (Mina et al. 2025) ---
LUNAR_RX_PARAMS = dict(
    T_eff=468.15, L_ad=0.0, L_pol=0.0, L_atm=0.0, Bn=2.0, Bp=12.0, T=0.01, b=10.0, D=0.1, Bf=0.2
)

# --- reference-sim transmit config: azimuth-averaged ACE block patterns + lower
#     per-block power [dBW]. The native mode instead keeps LuPNT's per-SVN 3-D LM/ACE
#     patterns and default power (loaded by setup_transmitters). ---
_REF_GPS_PTX_DBW = {"IIR": 15.0, "IIR_M": 14.3, "IIF": 14.3, "III": 14.3, "IIA": 14.3}
_REF_GAL_PTX_DBW = 14.0
_GPS_BLOCK_ACE = {
    "IIR": "Block-IIR_ACE",
    "IIR_M": "Block-IIR-M_ACE",
    "IIF": "Block-IIF_ACE",
    "III": "Block-IIF_ACE",
    "IIA": "Block-IIF_ACE",
}

_CONST = {  # name -> (SP3 prefix, GnssConst, primary frequency)
    "GPS": ("G", pnt.GnssConst.GPS, pnt.GnssFreq.L1),
    "GALILEO": ("E", pnt.GnssConst.GALILEO, pnt.GnssFreq.E1),
}
_WEEK_S = 7 * 86400.0


def _gps_prn_to_block() -> dict:
    """PRN -> block name (IIR / IIR_M / IIF / III / IIA) from the LuPNT gps_table.csv."""
    path = pnt.get_file_path("gps_table.csv")
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            out[int(row["PRN"])] = row["blockName"].strip()
    return out


class LunarGnssConfig:
    """Configuration builder for a lunar GNSS measurement simulation.

    Parameters
    ----------
    mode : {"native", "reference"}
        Transmit model (see module docstring).
    constellations : sequence of {"GPS", "GALILEO"}
        Which GNSS systems to include.
    receiver_antenna : str
        Stored antenna name for the lunar receiver (default ``"moongpsr"``).
    prop_step_s : float
        Fixed RK8 step [s] for the full-dynamics orbit propagation.
    receiver_error_terms : dict or None
        Optional broadcast-ephemeris/clock and oscillator/vibration 1-sigma terms
        (keys ``sigma_pr_eph_m``, ``sigma_pr_clk_m``, ``allan_deviation``, ``sigma_vib_deg``)
        root-sum-squared into the reported measurement sigmas. Default ``None`` -> all 0
        (a "true" sim where truth == receiver knowledge). Populate these for a future-epoch
        run in which no SP3/broadcast pair realizes the ephemeris/clock error naturally.
    earth_occult_alt_km : float
        Earth-limb altitude mask [km]. A line of sight is dropped when it passes within this
        altitude of the Earth's surface, i.e. the Earth occulting body radius is set to
        ``R_EARTH + earth_occult_alt_km``. Default ``0`` masks only rays blocked by the solid
        Earth; raise it to also reject rays that graze the ionosphere/plasmasphere and would
        otherwise carry unmodeled first-order delay.
    """

    def __init__(
        self,
        mode="native",
        constellations=("GPS", "GALILEO"),
        receiver_antenna="moongpsr",
        prop_step_s=300.0,
        receiver_error_terms=None,
        earth_occult_alt_km=0.0,
    ):
        if mode not in ("native", "reference"):
            raise ValueError(f"mode must be 'native' or 'reference', got {mode!r}")
        bad = [c for c in constellations if c not in _CONST]
        if bad:
            raise ValueError(f"unknown constellation(s) {bad}; choose from {list(_CONST)}")
        self.mode = mode
        self.constellations = tuple(constellations)
        self.receiver_antenna = receiver_antenna
        self.prop_step_s = float(prop_step_s)
        self.receiver_error_terms = dict(receiver_error_terms or {})
        self.earth_occult_alt_km = float(earth_occult_alt_km)
        self.constellation_objects = {}  # name -> GnssConstellation (filled by build)

    # -- full-dynamics orbit propagator (Earth 8x8 + Sun + Moon) ------------------------
    def _make_dynamics(self):
        dyn = pnt.NBodyDynamics()
        dyn.set_frame(pnt.GCRF)
        dyn.set_units(pnt.SI_UNITS)
        dyn.add_body(pnt.create_body(pnt.EARTH, 8, 8))
        dyn.add_body(pnt.create_body(pnt.SUN))
        dyn.add_body(pnt.create_body(pnt.MOON))
        dyn.set_autodiff(False)
        dyn.set_integrator(pnt.IntegratorType.RK8)
        dyn.set_time_step(self.prop_step_s)
        return dyn

    @staticmethod
    def _latest_sp3_seed(before_tai):
        """Latest available precise SP3, stepping back one week at a time from ``before_tai``."""
        t = float(before_tai)
        for _ in range(60):  # up to ~1 year back
            try:
                path = pnt.Sp3Loader.download_file_for_epoch(t, pnt.Time.TAI)
                return t, pnt.Sp3Loader([path])
            except Exception:
                t -= _WEEK_S
        raise RuntimeError("no precise SP3 product found within a year before the target epoch")

    def _receiver_params(self):
        rp = pnt.GnssReceiverParams()
        for k, v in LUNAR_RX_PARAMS.items():
            setattr(rp, k, v)
        for k, v in self.receiver_error_terms.items():
            setattr(rp, k, v)
        return rp

    def _build_constellation(self, name, epochs_tai, seed_tai, sp3, dyn):
        """Seed each PRN from SP3 and propagate (full dynamics) onto the epoch grid."""
        prefix, gnss_enum, freq = _CONST[name]
        seed_tdb = float(pnt.convert_time(seed_tai, pnt.TAI, pnt.TDB))
        # Ephemeris grid = measurement epochs padded by a margin so the interpolant covers
        # the transmit times (receive - light-time, ~1.3 s Earth-Moon) at the grid edges.
        eph_tai = np.asarray(epochs_tai, float)
        pad = np.linspace(60.0, 1.0, 3)  # a few lead-in / lead-out knots
        eph_tai = np.unique(np.concatenate([eph_tai[0] - pad, eph_tai, eph_tai[-1] + pad]))
        grid_tdb = np.array([float(pnt.convert_time(float(t), pnt.TAI, pnt.TDB)) for t in eph_tai])
        t_arr = np.concatenate([[seed_tdb], grid_tdb])
        prns, rv_hist = [], []
        for prn in range(1, 37):
            sat = f"{prefix}{prn:02d}"
            if not sp3.has_satellite(sat):
                continue
            try:
                rv_ecef, _ = sp3.get_pos_vel_clock(sat, seed_tai)
            except RuntimeError:
                continue
            rv0 = np.array(
                pnt.convert_frame(seed_tdb, np.asarray(rv_ecef, float), pnt.ECEF, pnt.ECI)
            )
            traj = np.array(dyn.propagate(rv0, t_arr))  # [len(t_arr), 6] ECI
            prns.append(prn)
            rv_hist.append(traj[1:])  # grid states [len(eph_tai), 6]

        gc = pnt.GnssConstellation(gnss_enum)
        gc.set_satellite_states(prns, eph_tai, rv_hist)
        gc.setup_transmitters()  # native per-SVN patterns + power
        if self.mode == "reference":
            self._apply_reference_tx(gc, name, prns, freq)
        return gc, prns, freq

    def _apply_reference_tx(self, gc, name, prns, freq):
        """Swap to azimuth-averaged ACE block patterns + lower per-block power."""
        if name == "GPS":
            prn2blk = _gps_prn_to_block()
            for prn in prns:
                blk = prn2blk.get(prn)
                if blk is None:
                    continue
                gc.set_transmitter_antenna(prn, freq, pnt.Antenna(_GPS_BLOCK_ACE[blk]))
                gc.set_transmit_power_dbw(prn, freq, _REF_GPS_PTX_DBW[blk])
        else:  # Galileo already uses the Galileo_E1 pattern; only the power is lower
            for prn in prns:
                gc.set_transmit_power_dbw(prn, freq, _REF_GAL_PTX_DBW)

    def _measurement_options(self):
        opts = pnt.GnssMeasurementOptions()
        opts.frame = pnt.ECI
        opts.receive_time_scale = pnt.TAI
        opts.ephemeris_time_scale = pnt.TAI
        opts.solve_light_time = True
        opts.apply_visibility = True
        opts.apply_cn0_threshold = False
        opts.apply_ionosphere_plasma_delay = False
        opts.tx_gain_model = (
            pnt.GnssMeasurementOptions.TxGainModel.AZIMUTH_AVERAGED
            if self.mode == "reference"
            else pnt.GnssMeasurementOptions.TxGainModel.FULL_2D
        )
        return opts

    def build_measurements(self, epochs_tai, seed_tai=None):
        """Build a configured :class:`GNSSMeasurements` for the given epoch grid (TAI seconds).

        ``seed_tai`` is the SP3 seed epoch; if ``None``, the latest precise SP3 before the
        first epoch is used (stepping back weekly), so future epochs are handled.
        """
        epochs_tai = np.asarray(epochs_tai, float)
        if seed_tai is None:
            seed_tai, sp3 = self._latest_sp3_seed(epochs_tai[0])
        else:
            sp3 = pnt.Sp3Loader([pnt.Sp3Loader.download_file_for_epoch(seed_tai, pnt.Time.TAI)])
        self.seed_tai = float(seed_tai)

        dyn = self._make_dynamics()
        built = [
            (name, *self._build_constellation(name, epochs_tai, seed_tai, sp3, dyn))
            for name in self.constellations
        ]
        self.constellation_objects = {name: gc for name, gc, _, _ in built}

        (name0, gc0, _, freq0) = built[0]
        gm = pnt.GNSSMeasurements(gc0)
        gm.set_frequency(freq0)
        for name, gc, _, freq in built[1:]:
            gm.add_constellation(gc, freq)

        gm.set_options(self._measurement_options())
        gm.set_receiver_params(self._receiver_params())
        gm.set_receiver_antenna(pnt.Antenna(self.receiver_antenna))

        # Standard lunar setup: Earth + Moon occlude the line of sight; the receiver
        # boresight points at Earth (ECI origin); the Sun drives the transmitter yaw model.
        earth = pnt.GnssOccludingBody()
        earth.radius_m = pnt.R_EARTH + self.earth_occult_alt_km * 1e3
        moon = pnt.GnssOccludingBody()
        moon.radius_m = pnt.R_MOON
        moon.position_provider = lambda t: np.array(
            pnt.get_body_pos(
                float(pnt.convert_time(t, pnt.TAI, pnt.TDB)), pnt.EARTH, pnt.MOON, pnt.ECI
            )
        )
        gm.set_occluding_bodies([earth, moon])
        gm.set_sun_position_provider(
            lambda t: np.array(
                pnt.get_body_pos(
                    float(pnt.convert_time(t, pnt.TAI, pnt.TDB)), pnt.EARTH, pnt.SUN, pnt.ECI
                )
            )
        )
        gm.set_boresight_target_provider(lambda t: np.zeros(3))
        return gm


def lunar_gnss_measurements(epochs_tai, mode="native", seed_tai=None, **kwargs):
    """Convenience wrapper: build a lunar :class:`GNSSMeasurements` in one call.

    Equivalent to ``LunarGnssConfig(mode=mode, **kwargs).build_measurements(epochs_tai, seed_tai)``.
    """
    return LunarGnssConfig(mode=mode, **kwargs).build_measurements(epochs_tai, seed_tai)
