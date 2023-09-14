import pylupnt as pnt
import numpy as np
import sys
import pytest

sys.path.append("pylupntutil")
import pntautodiff as ad


def test_KeplerianDynamics_ClassicalOE():
    # Classical Orbital elements
    a = 6541.4
    e = 0.6
    i = np.deg2rad(65.5)
    Omega = np.deg2rad(90.0)
    w = np.deg2rad(0.0)
    M = np.deg2rad(0.0)

    mu = pnt.MU_MOON
    coe_state = pnt.ClassicalOE(a, e, i, Omega, w, M, cs=pnt.CoordSystem.MI)
    coe_analytical = coe_state.get_vector().asarray()

    # Keplerian dynamics
    kep_dyn = pnt.KeplerianDynamics(mu)

    # Test propagation
    dt = 10.0
    n = np.sqrt(mu / np.power(a, 3))
    for _ in range(100):
        kep_dyn.propagate(coe_state, dt)
        coe_analytical[5] = pnt.wrapToPi(coe_analytical[5] + n * dt)
        np.testing.assert_almost_equal(coe_state.get_vector().asarray(), coe_analytical)
        dt += 2.0

    # Test propagation with STM
    def propagate_with_stm_num(coe, dt):
        def f_prop(coe_vec):
            coe_tmp = pnt.ClassicalOE(ad.realvec(coe_vec.asarray()))
            kep_dyn.propagate(coe_tmp, dt)
            return coe_tmp.get_vector()

        # Numerical STM
        coe_vec = ad.realvec(coe.get_vector().asarray())
        stm_num = ad.numerical_jacobian(f_prop, coe_vec, eps=1e-6)
        # Propagate
        kep_dyn.propagate(coe, dt)
        return stm_num

    dt = 10.0
    coe_numerical = pnt.ClassicalOE(coe_state.get_vector())
    stm_analytical = np.eye(6, dtype=float)
    stm_dyn = np.zeros((6, 6), dtype=float)
    for _ in range(100):
        stm_dyn = kep_dyn.propagate_with_stm(coe_state, dt)
        stm_analytical[5, 0] = -3.0 / 2.0 * n / a * dt
        stm_num = propagate_with_stm_num(coe_numerical, dt)

        np.testing.assert_almost_equal(stm_dyn, stm_analytical)
        np.testing.assert_almost_equal(stm_dyn, stm_num)
        dt += 2.0


def test_CartesianTwoBodyDynamics():
    # Classical Orbital elements
    a = 6541.4  # [km]
    e = 0.6  # [-]
    i = np.deg2rad(65.5)  # [rad]
    Omega = np.deg2rad(90.0)  # [rad]
    w = np.deg2rad(0.0)  # [rad]
    M = np.deg2rad(0.0)  # [rad]

    mu = pnt.MU_MOON
    coe_state = pnt.ClassicalOE(a, e, i, Omega, w, M, cs=pnt.CoordSystem.MI)
    cart_state = pnt.coe_to_cart(coe_state, mu)
    cart_vector = ad.realvec(cart_state.get_vector().asarray(), "X")

    # Two body dynamics
    kep_dyn = pnt.KeplerianDynamics(mu)
    tb_dyn = pnt.CartesianTwoBodyDynamics(mu)

    # Propagation
    dt = 10.0
    for _ in range(5):
        kep_dyn.propagate(coe_state, dt)
        tb_dyn.propagate(cart_state, 0, dt, 1.0)
        # tb_dyn.propagate(cart_vector, 0, dt, 1.0)

        cart_vector_kep = pnt.coe_to_cart(coe_state.get_vector(), mu).asarray()
        np.testing.assert_array_almost_equal(
            cart_vector_kep, cart_state.get_vector().asarray()
        )
        # np.testing.assert_array_almost_equal(cart_vector_kep, cart_vector)

    # Propagation with STM
    def propagate_with_stm_num(cart, t0, tf, dt):
        def f_prop(cart_vec):
            cart_tmp = pnt.CartesianState(ad.realvec(cart_vec.asarray()))
            tb_dyn.propagate(cart_tmp, t0, tf, dt)
            return cart_tmp.get_vector()

        # Numerical STM
        cart_vec = ad.realvec(cart.get_vector().asarray())
        stm_num = ad.numerical_jacobian(f_prop, cart_vec, eps=1e-5)
        # Propagate
        tb_dyn.propagate(cart, t0, tf, dt)
        return stm_num

    dt = 10.0
    cart_numerical = pnt.CartesianState(cart_state.get_vector())
    stm_dyn = np.zeros((6, 6), dtype=float)
    for _ in range(5):
        stm_dyn = tb_dyn.propagate_with_stm(cart_state, 0, dt, 1.0)
        stm_numerical = propagate_with_stm_num(cart_numerical, 0, dt, 1.0)

        np.testing.assert_array_almost_equal(
            cart_state.get_vector().asarray(), cart_numerical.get_vector().asarray()
        )
        np.testing.assert_array_almost_equal(stm_dyn, stm_numerical)


if __name__ == "__main__":
    test_KeplerianDynamics_ClassicalOE()
    test_CartesianTwoBodyDynamics()
    print("All tests passed!")
