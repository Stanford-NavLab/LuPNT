import pylupnt as pnt
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

# my modules
try:
    from .orbit_manager import OrbitManager
    from .ephemeris_sim import EphemerisSimulation
except ImportError:
    from orbit_manager import OrbitManager
    from ephemeris_sim import EphemerisSimulation


orbits = ["Polar"]  # ["Moonlight", "LNSS"] # ["ELFO", "Polar", "NRHO"]
poly_types = ["chebyshev"]
use_meq = False
stop_order_iter = False  # stop

overwrite_orbitsetup = False
overwrite_fit = True
overwrite_datasize = True

for orbit in orbits:
    print("---------------------------------------------------")
    print("Running ephemeris simulation for {}".format(orbit))
    print("---------------------------------------------------")

    # time settings
    if (
        orbit == "ELFO"
        or orbit == "Polar"
        or orbit == "CLFO"
        or orbit == "LCRNS"
        or orbit == "LNSS"
    ):
        dt = 0.1
        n_period = 3
        sphm = [80, 80]  # Spherical harmonic model degree and order
    else:
        dt = 10.0
        n_period = 2
        sphm = [8, 8]  # Spherical harmonic model degree and order

    if orbit == "Polar":
        max_order_default = 30
    else:
        max_order_default = 20

    add_earth = True
    add_sun = True

    pnt.set_lupnt_epoch(0.0)

    basedir = "/Users/keidaiiiyama/Documents/sw_navlab/LuPNT-private/projects/Ephemeris/"
    orbm_save_dir = basedir + "data/orbits"

    orbm = OrbitManager(
        orbit, dyn=None, n_period=n_period, dt=dt, data_dir=orbm_save_dir, overwrite=False
    )

    # setup ephemeris simulation
    esim_dir = basedir + "data/ephemeris/"
    fit_mins = [60, 120, 240, 360, 480]
    dt_fit = 60.0
    dt_eval = 1.0

    # setup ephemeris configuration
    for fit_min in fit_mins:
        print(" ")
        print("--------------------------------")
        print("[fit_min: {}]".format(fit_min))
        esim = EphemerisSimulation(data_dir=esim_dir)

        for use_cheby_sampling in [True, False]:
            print(" [use_cheby_sampling: {}]".format(use_cheby_sampling))
            esim.setup_orbit(
                orbm,
                sample_M=30,
                fit_mins=[fit_min],
                use_cheby_sampling=use_cheby_sampling,
                dt_fit=dt_fit,
                dt_eval=dt_eval,
                overwrite=overwrite_orbitsetup,
            )

            for use_kep in [True, False]:
                if use_kep:
                    use_fourier_list = [True, False]
                else:  # cannot use Fourier for non-keplarian ephemeris
                    use_fourier_list = [False]

                for poly_type in poly_types:

                    if poly_type == "monomial":
                        max_order = 4  # bigger than this will cause numerical instability
                    else:
                        max_order = max_order_default

                    for use_fourier in use_fourier_list:
                        cleared = False
                        cleared_prev = False
                        order = 2  # minimum order is 2

                        while not cleared:
                            config = {
                                "order": order,
                                "use_kep": use_kep,
                                "use_rsw": False,
                                "use_fourier": use_fourier,
                                "use_meq": use_meq,
                                "poly_type": poly_type,
                                "sampling_type": "cheby" if use_cheby_sampling else "uniform",
                            }

                            ephem_config = esim.fit_ephemeris(
                                ephem_type="cartesian",
                                config=config,
                                print_errors=False,
                                fit_obj="lsq-cvx",
                                print_opt_results=False,
                                overwrite=overwrite_fit,
                            )

                            esim.compute_datasize(
                                configs=[ephem_config],
                                fit_mins=[fit_min],
                                precision=1e-2,  # 1 cm, 0.01 mm/s
                                debug=False,
                                overwrite=overwrite_datasize,
                            )

                            p95_pos = esim.fit_results[ephem_config][fit_min]["pos_p95"]
                            p95_vel = esim.fit_results[ephem_config][fit_min]["vel_p95"]
                            databit = esim.datasizes[ephem_config][fit_min]["total_bits"]

                            print(
                                "    type: {0} | order: {1}, kep: {2}  fourier: {3} | p95_pos: {4:.5f} m, p95_vel: {5:.5f} mm/s | databit: {6} bits".format(
                                    poly_type,
                                    order,
                                    use_kep,
                                    use_fourier,
                                    p95_pos,
                                    p95_vel,
                                    databit,
                                )
                            )

                            if stop_order_iter:
                                if p95_pos <= 10 and p95_vel <= 2.5 and databit <= 900:
                                    if cleared_prev:  # already cleared in previous order -> break
                                        cleared = True
                                        print(
                                            "  ====== Cleared Twice: Order: {}, Fourier: {} ========".format(
                                                order, use_fourier
                                            )
                                        )
                                        break
                                    else:  # satisfies constraints -> try one more order
                                        cleared_prev = True
                                        order += 1
                                elif databit > 900:
                                    print(
                                        "  ====== Reached databit > 900: Terminate fitting ======="
                                    )
                                    cleared = False
                                    break
                                else:  # does not satisfy constraints -> try one more order
                                    cleared = False
                                    order += 1
                            else:
                                if order == max_order:
                                    cleared = True
                                    break
                                order += 1

                        print(" ")
