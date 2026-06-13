.. _python_basic_tutorial:

Python ELFO Propagation and Frame Conversion
============================================

This example mirrors the C++ tutorial with the Python bindings. It defines an
ELFO with classical orbital elements in the Moon orbital plane frame, converts
the initial state to Moon-centered inertial coordinates, propagates it with
third-body perturbations, and converts the propagated states to the Moon
principal-axes frame.

.. code-block:: python

    import numpy as np
    import pylupnt as pnt

    t0_tdb = pnt.convert_time(
        pnt.gregorian_to_time(2030, 1, 1, 12, 0, 0),
        pnt.UTC,
        pnt.TDB,
    )

    coe_op = np.array(
        [
            6541.4e3,       # semi-major axis [m]
            0.60,           # eccentricity [-]
            65.5 * pnt.RAD, # inclination [rad]
            0.0 * pnt.RAD,  # right ascension of ascending node [rad]
            90.0 * pnt.RAD, # argument of periapsis [rad]
            0.0 * pnt.RAD,  # true anomaly [rad]
        ]
    )

    rv0_op = pnt.classical_to_cart(coe_op, pnt.GM_MOON)
    rv0_ci = pnt.convert_frame(t0_tdb, rv0_op, pnt.MOON_OP, pnt.MOON_CI)

    dynamics = pnt.NBodyDynamics()
    dynamics.set_frame(pnt.MOON_CI)
    dynamics.add_body(pnt.Body.Moon(20, 20))
    dynamics.add_body(pnt.Body.Earth())
    dynamics.add_body(pnt.Body.Sun())

    tspan = np.arange(t0_tdb, t0_tdb + 7.0 * pnt.SECS_DAY, pnt.SECS_HOUR)
    rv_ci = dynamics.propagate(rv0_ci, tspan)
    rv_pa = pnt.convert_frame(tspan, rv_ci, pnt.MOON_CI, pnt.MOON_PA)

    print("Initial Moon-CI state [m, m/s]:", rv0_ci)
    print("Final Moon-PA state [m, m/s]:", rv_pa[-1])

The dynamics model is configured in Moon-CI because the numerical propagator
expects a single integration frame. The frame conversion calls handle the
initial geometry and any output products that are more useful in a body-fixed
or orbit-related frame.
