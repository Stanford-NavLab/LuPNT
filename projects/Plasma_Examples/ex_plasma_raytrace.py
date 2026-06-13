import numpy as np
from pylupnt import plasma as tec
import sys

# visibility parameters
mainlobe_angle = 20.0  # mainlobe angle in degrees (only mainlobe is used)
min_h = 300  # minimum height of the receiver in km
max_min_alt = 1500  # maximum minimum altitude of the receiver in km
num_Omega = 72  # number of Omega satellites

# other parameters
use_moon = True  # if True, use moon as a receiver, if False use GEO satlellite as receiver
debug_prop = False  # if True, print propagation debug info
debug_corr = False  # correction disabled; no corrector output to print

# GNSS constellation
gps_sats = tec.setup_gnss_constellation("gps_2025_01_01.txt")
epoch_utc = gps_sats[0].epoch_utc_
num_gps = len(gps_sats)

# GEO satellite
id_geo = 999
a_geo = 42164.0
e_geo = 0.0001
inc_geo = 0.0001

# Moon
a_moon = 384400.0
i_moon = 23.44 * tec.DEG2RAD
e_moon = 0.0549

# IRI model
tec.set_iri_model("IRI2007")

# Store Tx and Rx conditions where ray path ionosphere
epoch_utcs = []
pos_txs = []
pos_rxs = []
min_alts = []
prns = []
rv_sats = []

case_idx = 0

for omi in range(num_Omega):
    Omega = tec.DEG2RAD * omi * (360 / num_Omega)

    if use_moon:
        a_rx = a_moon
        e_rx = e_moon
        inc_rx = i_moon
    else:
        a_rx = a_geo
        e_rx = e_geo
        inc_rx = inc_geo

    rx_sat = tec.Satellite(
        id_geo, np.array([a_rx, e_rx, inc_rx, Omega, 0.0, 0.0]), epoch_utc, tec.GM_EARTH
    )

    for i, sat in enumerate(gps_sats):
        pos_gps = sat.get_pos()
        pos_rx = rx_sat.get_pos()

        vis = tec.compute_vis(pos_gps, pos_rx, tec.RE + min_h, mainlobe_angle)
        min_alt = tec.compute_min_altitude(pos_gps, pos_rx, tec.RE)

        if vis and min_alt <= max_min_alt:
            print(
                f"Case: {case_idx} | Omega: {omi}/{num_Omega}, Omega: {Omega * 180 / np.pi:.2f} deg  PRN: {sat.id_} | Visibility: {'Yes' if vis else 'No'}, Minimum Altitude: {min_alt:.2f} km"
            )
            pos_tx = tec.solve_lt(sat, pos_rx, epoch_utc)

            epoch_utcs.append(epoch_utc)
            pos_txs.append(pos_tx)
            pos_rxs.append(pos_rx)
            min_alts.append(min_alt)
            prns.append(sat.id_)
            rv_sats.append(sat.posvel_)

            case_idx += 1

# Ray tracing simulation settings
# Ray tracing configuration
config = tec.RayTraceConfig()
config.freq_Hz = tec.freq_L1
config.step_size = 50.0  # step size in km
config.correction = True
config.fine_correction = True
config.cutoff_r = 4 * tec.RE  # > GPS orbit (~26,560 km); 5·RE ≈ 31,891 km
config.gradn_dx = 1.0  # gradient step size in km
config.integ_method = "RK4"  # integration method
config.kp = -1  # Kp index, -1 means it will be computed automatically
config.correction_method = "neldermead"
config.use_fortran_gcpm = True
config.corr_tol = 1.0  # correction tolerance in meters
config.use_adaptive_step = True

# Get Kp index if needed
datetime = tec.mjd_to_datetime(tec.tj2000_to_mjd(epoch_utc))
if config.kp < 0:
    config.kp = tec.get_kp_index(datetime)

print(f"Kp index: {config.kp}")
print(
    f"Date: Year: {datetime.year}, DOY: {datetime.doy}, Hour: {datetime.hour}, Minute: {datetime.min}, Second: {datetime.sec}"
)


def run_raytrace(idx):
    epoch_utc = epoch_utcs[idx]
    pos_tx = pos_txs[idx]
    pos_rx = pos_rxs[idx]

    sys.stdout.flush()

    print(f"Solving ray trace for PRN {gps_sats[0].id_} at epoch {epoch_utc} UTC")
    print(f"Tx position: {pos_tx}, Rx position: {pos_rx}")

    # Perform ray tracing
    pp = tec.trace_ray(epoch_utc, pos_tx, pos_rx, config, debug_prop=False, debug_corr=True)

    # Print results
    sys.stdout.flush()

    print(" ")
    print(f"[Raytrace Result] PRN = {i}")
    print(f"  Minimum Altitude: {min_alts[idx]} km")
    print(f"  TECU: {pp.tecu} TECU")
    print(f"  Total Delay: {pp.total_delay_m} m")
    # print(f"  Dist Total   : {pp.sf} m")
    # print(f"  Dist Straight: {pp.dist_straight_km} km")
    print(f"  Bend Delay    : {pp.dist_bend_m} m")
    print(f"  TEC  Delay    : {pp.tec_delay_m} m")
    print(f"  Final Pos Error: {np.linalg.norm(pp.corr_final_pos_err) * 1000} m")
    print(f"  Final Time Error: {pp.corr_final_time_err} s")

    return pp


# Case 1
print(" ")
print(" ")
max_alt_idx = np.argmax(min_alts)
print(
    f"Running ray trace for case with minimum altitude {min_alts[max_alt_idx]:.2f} km at index {max_alt_idx}"
)
pp1 = run_raytrace(max_alt_idx)
