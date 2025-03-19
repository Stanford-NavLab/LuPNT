import pylupnt as pnt
import numpy as np
import pytest


ABS_TOL = 1e-6
REL_TOL = 1e-6


def get_classical_oe():
    a = 5740  # [km] Semi-major axis
    e = 0.58  # [-] Eccentricity
    i = 54.9 * pnt.RAD  # [rad] Inclination
    O = 0.00 * pnt.RAD  # [rad] Right Ascension of Ascending Node
    w = 86.3 * pnt.RAD  # [rad] Argument of Perigee
    M = 0.00 * pnt.RAD  # [rad] True Anomaly

    coe = np.array([a, e, i, O, w, M])
    return coe


def get_classical_oe_mat(n):
    coe = get_classical_oe()
    coe_mat = np.tile(coe, (n, 1))
    coe_mat[:, 5] = np.linspace(0, 2 * np.pi, n)
    return coe_mat


def test_twobody_dynamics():
    # Constants
    N_sat = 3  # Number of satellites
    J2 = 0  # [-] J2 coefficient
    GM = pnt.GM_MOON  # [km^3/s^2] Gravitational parameter
    R_body = pnt.R_MOON  # [km] Radius of the central body

    # Time
    N_steps = 4
    dt = 10  # [s] Integration time step
    Dt = 6 * pnt.SECS_HOUR  # [s] Total integration time
    t0 = pnt.gregorian2time(2024, 6, 1, 12, 45, 30)  # [s] Start time, TAI
    tspan = np.linspace(0, Dt, N_steps)  # [s] Time span
    tfs = t0 + tspan  # [s] Times, TAI

    # Dynamics
    dyn_kep = pnt.KeplerianDynamics(GM)
    dyn_cart = pnt.CartesianTwoBodyDynamics(GM, pnt.IntegratorType.RK4)
    dyn_cart_j2 = pnt.J2CartTwoBodyDynamics(GM, J2, R_body, pnt.IntegratorType.RK4)
    dyn_kep_j2 = pnt.J2KeplerianDynamics(GM, J2, R_body, pnt.IntegratorType.RK4)

    # Time step
    dyn_cart.set_time_step(dt)
    dyn_cart_j2.set_time_step(dt)
    dyn_kep_j2.set_time_step(dt)

    # Vec6
    coe_vec = get_classical_oe()
    rv_vec = pnt.classical2cart(coe_vec, GM)
    coe_j2_vec = coe_vec
    rv_j2_vec = rv_vec

    # States
    coe_state = pnt.ClassicalOE(coe_vec)
    coe_j2_state = pnt.ClassicalOE(coe_j2_vec)
    rv_state = pnt.CartesianOrbitState(rv_vec)
    rv_j2_state = pnt.CartesianOrbitState(rv_j2_vec)

    # ************************************************************************************************
    # Multple times
    # ************************************************************************************************

    # (N_steps x 6)
    for shape in ("array", "column", "row"):
        if shape == "array":
            coe_prop = dyn_kep.propagate(coe_vec, t0, tfs)
            coe_j2_prop = dyn_kep_j2.propagate(coe_j2_vec, t0, tfs)
            rv_prop = dyn_cart.propagate(rv_vec, t0, tfs)
            rv_j2_prop = dyn_cart_j2.propagate(rv_j2_vec, t0, tfs)
        elif shape == "column":
            coe_prop = dyn_kep.propagate(coe_vec.reshape(6, 1), t0, tfs)
            coe_j2_prop = dyn_kep_j2.propagate(coe_j2_vec.reshape(6, 1), t0, tfs)
            rv_prop = dyn_cart.propagate(rv_vec.reshape(6, 1), t0, tfs)
            rv_j2_prop = dyn_cart_j2.propagate(rv_j2_vec.reshape(6, 1), t0, tfs)
        elif shape == "row":
            coe_prop = dyn_kep.propagate(coe_vec.reshape(1, 6), t0, tfs)
            coe_j2_prop = dyn_kep_j2.propagate(coe_j2_vec.reshape(1, 6), t0, tfs)
            rv_prop = dyn_cart.propagate(rv_vec.reshape(1, 6), t0, tfs)
            rv_j2_prop = dyn_cart_j2.propagate(rv_j2_vec.reshape(1, 6), t0, tfs)

    coe_j2_prop[:, 5] = pnt.wrap2pi(coe_j2_prop[:, 5])

    # Check rows
    assert coe_prop.shape[0] == N_steps
    assert coe_j2_prop.shape[0] == N_steps
    assert rv_prop.shape[0] == N_steps
    assert rv_j2_prop.shape[0] == N_steps

    # Check values
    np.testing.assert_allclose(coe_prop, coe_j2_prop, atol=ABS_TOL)
    np.testing.assert_allclose(coe_prop, pnt.cart2classical(rv_prop, GM), atol=ABS_TOL)
    np.testing.assert_allclose(
        coe_prop, pnt.cart2classical(rv_j2_prop, GM), atol=ABS_TOL
    )

    # (N_steps x 6)
    coe_prop = dyn_kep.propagate(coe_vec, t0, tfs)
    coe_j2_prop = dyn_kep_j2.propagate(coe_j2_vec, t0, tfs)
    rv_prop = dyn_cart.propagate(rv_vec, t0, tfs)
    rv_j2_prop = dyn_cart_j2.propagate(rv_j2_vec, t0, tfs)

    coe_j2_prop[:, 5] = pnt.wrap2pi(coe_j2_prop[:, 5])

    # Check rows
    assert coe_prop.shape[0] == N_steps
    assert coe_j2_prop.shape[0] == N_steps
    assert rv_prop.shape[0] == N_steps
    assert rv_j2_prop.shape[0] == N_steps

    # Check values
    np.testing.assert_allclose(coe_prop, coe_j2_prop, atol=ABS_TOL)
    np.testing.assert_allclose(coe_prop, pnt.cart2classical(rv_prop, GM), atol=ABS_TOL)
    np.testing.assert_allclose(
        coe_prop, pnt.cart2classical(rv_j2_prop, GM), atol=ABS_TOL
    )

    # ************************************************************************************************
    # Multple vectors
    # ************************************************************************************************

    # (N_sat x 6)
    coe_prop_sat = get_classical_oe_mat(N_sat)
    coe_j2_prop_sat = coe_prop_sat
    rv_matx6 = pnt.classical2cart(coe_prop_sat, GM)
    rv_j2_prop_sat = pnt.classical2cart(coe_prop_sat, GM)

    # Check times
    np.testing.assert_allclose(t0, tfs[0], atol=ABS_TOL)

    # Propagate loop
    for i in range(1, N_steps):
        coe_prop_sat = dyn_kep.propagate(coe_prop_sat, tfs[i - 1], tfs[i])
        coe_j2_prop_sat = dyn_kep_j2.propagate(coe_j2_prop_sat, tfs[i - 1], tfs[i])
        rv_matx6 = dyn_cart.propagate(rv_matx6, tfs[i - 1], tfs[i])
        rv_j2_prop_sat = dyn_cart_j2.propagate(rv_j2_prop_sat, tfs[i - 1], tfs[i])

        coe_j2_prop_sat[:, 5] = pnt.wrap2pi(coe_j2_prop_sat[:, 5])

        # Check values
        np.testing.assert_allclose(coe_prop_sat[0], coe_prop[i], atol=ABS_TOL)
        np.testing.assert_allclose(coe_prop_sat, coe_j2_prop_sat, atol=ABS_TOL)
        np.testing.assert_allclose(
            coe_prop_sat, pnt.cart2classical(rv_matx6, GM), atol=ABS_TOL
        )
        np.testing.assert_allclose(
            coe_prop_sat, pnt.cart2classical(rv_j2_prop_sat, GM), atol=ABS_TOL
        )

    # ************************************************************************************************
    # Vectors
    # ************************************************************************************************

    for shape in ("array", "column", "row"):
        coe_vec_prop = coe_vec
        coe_j2_vec_prop = coe_j2_vec
        rv_vec_prop = rv_vec
        rv_j2_vec_prop = rv_j2_vec

        for i in range(1, N_steps):
            if shape == "array":
                coe_vec_prop = dyn_kep.propagate(coe_vec_prop, tfs[i - 1], tfs[i])
                coe_j2_vec_prop = dyn_kep_j2.propagate(
                    coe_j2_vec_prop, tfs[i - 1], tfs[i]
                )
                rv_vec_prop = dyn_cart.propagate(rv_vec_prop, tfs[i - 1], tfs[i])
                rv_j2_vec_prop = dyn_cart_j2.propagate(
                    rv_j2_vec_prop, tfs[i - 1], tfs[i]
                )
            elif shape == "column":
                coe_vec_prop = dyn_kep.propagate(
                    coe_vec_prop.reshape(6, 1), tfs[i - 1], tfs[i]
                )
                coe_j2_vec_prop = dyn_kep_j2.propagate(
                    coe_j2_vec_prop.reshape(6, 1), tfs[i - 1], tfs[i]
                )
                rv_vec_prop = dyn_cart.propagate(
                    rv_vec_prop.reshape(6, 1), tfs[i - 1], tfs[i]
                )
                rv_j2_vec_prop = dyn_cart_j2.propagate(
                    rv_j2_vec_prop.reshape(6, 1), tfs[i - 1], tfs[i]
                )
            elif shape == "row":
                coe_vec_prop = dyn_kep.propagate(
                    coe_vec_prop.reshape(1, 6), tfs[i - 1], tfs[i]
                )
                coe_j2_vec_prop = dyn_kep_j2.propagate(
                    coe_j2_vec_prop.reshape(1, 6), tfs[i - 1], tfs[i]
                )
                rv_vec_prop = dyn_cart.propagate(
                    rv_vec_prop.reshape(1, 6), tfs[i - 1], tfs[i]
                )
                rv_j2_vec_prop = dyn_cart_j2.propagate(
                    rv_j2_vec_prop.reshape(1, 6), tfs[i - 1], tfs[i]
                )

            # coe_state = pnt.KeplerianDynamics.propagate(coe_state, tfs[i - 1], tfs[i])
            # coe_j2_state = pnt.J2KeplerianDynamics.propagate(
            #     coe_j2_state, tfs[i - 1], tfs[i]
            # )
            # rv_state = pnt.CartesianTwoBodyDynamics.propagate(
            #     rv_state, tfs[i - 1], tfs[i]
            # )
            # rv_j2_state = pnt.J2CartTwoBodyDynamics.propagate(
            #     rv_j2_state, tfs[i - 1], tfs[i]
            # )

            coe_j2_vec_prop[5] = pnt.wrap2pi(coe_j2_vec_prop[5])

            # Check values
            np.testing.assert_allclose(coe_vec_prop, coe_prop[i], atol=ABS_TOL)
            np.testing.assert_allclose(coe_vec_prop, coe_j2_vec_prop, atol=ABS_TOL)
            np.testing.assert_allclose(
                coe_vec_prop, pnt.cart2classical(rv_vec_prop, GM), atol=ABS_TOL
            )
            np.testing.assert_allclose(
                coe_vec_prop, pnt.cart2classical(rv_j2_vec_prop, GM), atol=ABS_TOL
            )


if __name__ == "__main__":
    np.set_printoptions(precision=3, suppress=True)
    test_twobody_dynamics()
