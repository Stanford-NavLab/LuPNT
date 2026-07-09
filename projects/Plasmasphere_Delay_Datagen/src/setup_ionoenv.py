import pylupnt as pnt
import numpy as np
import os
import tecsimpy as tec
import plotly.graph_objects as go

from src.gnss_meas import GNSSMeas


def setup_lcrns_sats(n_orbit=6.0, dt=1, savefig=False, overwrite=True):
    """
    Setup the LCRNS satellites.
    """
    # reference time for orbit definition and propagation
    refyear = 2027
    refmonth = 3
    refday = 1
    refhour = 0
    refminute = 0
    refsecond = 0

    # simulation time used for gnss measurement generation
    year = 2025
    month = 3
    day = 1
    hour = 12
    minute = 0
    second = 0

    # t0_tai = pnt.convert_time(pnt.gregorian2time(year, month, day, hour, minute, second), pnt.UTC, pnt.TAI)
    t0_tai_ref = pnt.convert_time(
        pnt.gregorian_to_time(refyear, refmonth, refday, refhour, refminute, refsecond),
        pnt.UTC,
        pnt.TAI,
    )
    t0_tai = pnt.convert_time(
        pnt.gregorian_to_time(year, month, day, hour, minute, second), pnt.UTC, pnt.TAI
    )

    # a, e, i, Omega, w, M0 in PA frame
    N_sc = 5  # number of spacecraft
    svoe = np.zeros((N_sc, 6))
    svoe[0] = np.array([11315.936501, 0.691982, 59.373229, 321.019197, 92.494031, 0.000000])
    svoe[1] = np.array([11317.948675, 0.691982, 58.951732, 320.997768, 92.505016, 180.000000])
    svoe[2] = np.array([11305.413654, 0.691982, 52.733096, 81.148790, 92.062891, 140.049207])
    svoe[3] = np.array([11326.302154, 0.691982, 52.513419, 81.138818, 92.068945, 195.992393])
    svoe[4] = np.array([11307.882863, 0.691982, 56.310396, 204.889626, 85.444071, 164.007607])
    svoe[:, 0] = svoe[:, 0] * 1000  # sma
    svoe[:, 2:] = np.deg2rad(svoe[:, 2:])

    print(f"t0 TAI: {t0_tai}")

    a = svoe[3, 0]  # semi-major axis (first one)
    t_orbit = 2 * np.pi * np.sqrt(a**3 / pnt.GM_MOON)  # orbital period
    tf = n_orbit * t_orbit
    N_t = int(tf / dt) + 1  # [-] Number of time steps
    tf = (N_t - 1) * dt  # [s] Final time adjusted to integer number of steps
    tspan = np.linspace(0, tf, N_t)  # [s] Time since first epoch
    t_tai = t0_tai + tspan  # [s] Time in TAI

    dt = tspan[1] - tspan[0]  # adjust dt
    print(f"Simulating for {tf/3600:.2f} hours ({N_t} steps of {dt:.6f} sec)")

    # Dynamics ----------------------------------------------------------------
    pnt.set_lupnt_epoch(0)
    dyn = pnt.NBodyDynamics()
    dyn.set_integrator(pnt.IntegratorType.RKF45)
    dyn.set_integrator_params(pnt.IntegratorParams(max_iter=100, abstol=1e-12, reltol=1e-12))
    dyn.add_body(pnt.Body.Moon(50, 50))
    dyn.add_body(pnt.Body.Earth())
    dyn.add_body(pnt.Body.Sun())
    dyn.set_frame(pnt.MOON_CI)
    dyn.set_srp_coeff(CR=1.8, area=1.0, mass=850.0)
    dyn.set_time_step(dt)

    print("SRP Coeff: ", dyn.get_srp_coeff())

    # convert to rv and propagate
    rv_m2sc_ci = np.zeros((N_sc, N_t, 6))
    rv_e2sc_ecef = np.zeros((N_sc, N_t, 6))

    basepath = pnt.get_output_dir()
    orbdir = os.path.join(basepath, "iono_delay", "orbits")
    if not os.path.exists(orbdir):
        os.makedirs(orbdir)

    lcrns_file = orbdir + "/lcrns_{0:d}_{1:d}_{2:d}_{3:d}_norbit_{4:.0f}_dt_{5:.1f}.npz".format(
        year, month, day, hour, n_orbit, dt
    )

    if os.path.exists(lcrns_file) and not overwrite:
        print(f"Loading existing data from {lcrns_file}...")
        data = np.load(lcrns_file)
        t_tai = data["t_tai"]
        rv_m2sc_ci = data["rv_m2sc_ci"]
        rv_e2sc_ecef = data["rv_e2sc_ecef"]

        rv0_mci = rv_m2sc_ci[:, 0, :]
        for i in range(N_sc):
            rv0_op = pnt.convert_frame(
                t0_tai_ref, rv0_mci[i], pnt.MOON_CI, pnt.MOON_OP, rotate_only=True
            )
            coe_op = pnt.cart_to_classical(rv0_op, pnt.GM_MOON)
            print(f"Initial orbital elements (LUPNT frame) for spacecraft {i}:")
            print(" a (km): ", coe_op[0] / 1000)
            print(" e      : ", coe_op[1])
            print(" i (deg): ", np.rad2deg(coe_op[2]))
            print(" Omega (deg): ", np.rad2deg(coe_op[3]))
            print(" w (deg): ", np.rad2deg(coe_op[4]))
            print(" M0 (deg): ", np.rad2deg(coe_op[5]))
            print(" ")
    else:
        print(f"File {lcrns_file} not found. Propagating orbits...")

        for i in range(N_sc):
            print(f"Propagating orbit {i + 1}/{N_sc}...")
            # convert to classical orbital elements
            rv_tref_pa = pnt.classical_to_cart(svoe[i], pnt.GM_MOON)
            rv_tref_op = pnt.convert_frame(
                t0_tai_ref, rv_tref_pa, pnt.MOON_PA, pnt.MOON_OP, rotate_only=True
            )

            coe_op = pnt.cart_to_classical(rv_tref_op, pnt.GM_MOON)

            print(f"Initial orbital elements (LUPNT frame) for spacecraft {i}:")
            print(" a (km): ", coe_op[0] / 1000)
            print(" e      : ", coe_op[1])
            print(" i (deg): ", np.rad2deg(coe_op[2]))
            print(" Omega (deg): ", np.rad2deg(coe_op[3]))
            print(" w (deg): ", np.rad2deg(coe_op[4]))
            print(" M0 (deg): ", np.rad2deg(coe_op[5]))
            print(" ")

            # Use the same Op frame state for the propagation
            rv_t0_ci = pnt.convert_frame(
                t0_tai, rv_tref_op, pnt.MOON_OP, pnt.MOON_CI, rotate_only=True
            )  # convert to LUPNT frame

            # propagate the orbit
            rv_m2sc_ci[i] = dyn.propagate(rv_t0_ci, t_tai)

        # convert to ECEF frame
        for i in range(N_sc):
            rv_e2sc_ecef[i] = pnt.convert_frame(
                t_tai, rv_m2sc_ci[i], pnt.MOON_CI, pnt.ECEF, rotate_only=False
            )

        # add the south pole
        rv_south_pole_pa = np.zeros((N_t, 6))
        rv_south_pole_pa[:, 2] = -pnt.R_MOON
        rv_south_pole_ci = pnt.convert_frame(
            t_tai, rv_south_pole_pa, pnt.MOON_PA, pnt.MOON_CI
        )  # [N, 6]
        print("South Pole position (CI) at t0:", rv_south_pole_ci[0])

        # add to the rv_m2sc_ci
        rv_m2sc_ci = np.concatenate((rv_m2sc_ci, rv_south_pole_ci[np.newaxis, :, :]), axis=0)
        N_sc += 1  # add south pole
        rv_e2sc_ecef = np.concatenate(
            (
                rv_e2sc_ecef,
                pnt.convert_frame(
                    t_tai, rv_south_pole_ci, pnt.MOON_CI, pnt.ECEF, rotate_only=False
                )[np.newaxis, :, :],
            ),
            axis=0,
        )

        # save the propagated orbits to file
        np.savez(lcrns_file, t_tai=t_tai, rv_m2sc_ci=rv_m2sc_ci, rv_e2sc_ecef=rv_e2sc_ecef)

    if savefig:
        fig = go.Figure()
        pnt.plot.plot_orbits(fig, rv_m2sc_ci)  # [N, t, 3]
        pnt.plot.plot_body(
            fig,
            pnt.MOON,
            size_factor=2,
            alpha=0.5,
        )
        pnt.plot.set_view(fig, -80, 20, 2.5)
        fig.update_layout(showlegend=True, width=400, height=400)

        # save the figure as pdf
        if not os.path.exists(orbdir + "/figures"):
            os.makedirs(orbdir + "/figures")

        fig.write_image(orbdir + "/figures/lcrns_orbits.pdf")
        fig.show()

    return t_tai, rv_m2sc_ci, rv_e2sc_ecef, N_sc


def setup_gnss_constellation(
    t_tai,
    rv_m2sc_ci,
    rv_e2sc_ecef=None,
    gps_datetime=None,
    sp3_prop_method="interp",
    savefig=False,
    overwrite_orbit=False,
    overwrite_measurements=False,
    consider_faults=False,
):
    """
    Setup the GNSS satellites.
    """
    gnss_meas_dir = os.path.join(pnt.get_output_dir(), "iono_delay")
    gnss_meas = GNSSMeas(
        t_tai,
        rv_m2sc_ci,
        basepath=gnss_meas_dir,
        consider_faults=consider_faults,
        rv_e2sc_ecef=rv_e2sc_ecef,
    )
    gnss_meas.setup_gnss(
        gnss_consts=["GPS", "GALILEO", "QZSS"],
        gps_datetime=gps_datetime,
        overwrite=overwrite_orbit,
        sp3_prop_method=sp3_prop_method,
    )
    gnss_meas.setup_measurements(cn0_threshold=15.0, overwrite=overwrite_measurements)
    sat_labels = ["LCRNS_1", "LCRNS_2", "LCRNS_3", "LCRNS_4", "LCRNS_5", "South Pole"]

    if savefig:
        if not os.path.exists(gnss_meas_dir + "/figures"):
            os.makedirs(gnss_meas_dir + "/figures")

        gnss_figname = os.path.join(gnss_meas_dir, "figures/gnss_sats.pdf")
        fig = gnss_meas.plot_gnss_orbit(savefig=True, filename=gnss_figname)

        gnss_num_filename = os.path.join(gnss_meas_dir, "figures/gnss_num_tracked_sats.pdf")
        fig2 = gnss_meas.plot_num_tracked_sats(sat_labels, savefig=True, filename=gnss_num_filename)

    return gnss_meas


def dynamics_with_stm(x):

    x = x[:n]  # state vector
    stm_flat = x[n:]  # flattened STM
    stm = stm_flat.reshape((n, n))  # reshape to n x n

    # state derivative
    dxdt = dynamics(x)

    # STM derivative
    A = compute_jacobian(x)  # compute the Jacobian of the dynamics

    dstm_dt = A @ stm  # STM derivative

    return np.concatenate((dxdt, dstm_dt.flatten()))  # concatenate state and STM derivatives
