"""Unit tests for the `pylupnt` orbit-dynamics bindings:
`KeplerianDynamics`, `CartesianTwoBodyDynamics`, `JToCartTwoBodyDynamics`,
`J2KeplerianDynamics`.

Propagation over a time grid uses the ``propagate(x0, tfs)`` overload, where
``tfs`` is the vector of absolute epochs (its first element is the start).
"""

import numpy as np
import pylupnt as pnt
import pytest

ABS_TOL = 1e-6
GM = pnt.GM_MOON
R_BODY = pnt.R_MOON

# A moderately eccentric lunar orbit: [a, e, i, RAAN, argp, M] (SI + radians).
COE0 = np.array([5740e3, 0.3, np.radians(54.9), 0.0, np.radians(86.3), 0.0])

T0 = pnt.convert_time(pnt.gregorian_to_time(2026, 1, 1, 0, 0, 0), pnt.Time.TDB, pnt.Time.TAI)
TFS = T0 + np.linspace(0.0, 3 * pnt.SECS_HOUR, 5)


def _specific_energy(rv):
    return 0.5 * np.dot(rv[3:], rv[3:]) - GM / np.linalg.norm(rv[:3])


# ---------------------------------------------------------------------------
# Keplerian vs Cartesian two-body agreement
# ---------------------------------------------------------------------------


def test_keplerian_matches_cartesian_two_body():
    rv0 = np.asarray(pnt.classical_to_cart(COE0, GM))

    kep = pnt.KeplerianDynamics(GM)
    cart = pnt.CartesianTwoBodyDynamics(GM)
    cart.set_time_step(10.0)

    coe_prop = np.asarray(kep.propagate(COE0, TFS))
    rv_prop = np.asarray(cart.propagate(rv0, TFS))
    assert coe_prop.shape == (len(TFS), 6)
    assert rv_prop.shape == (len(TFS), 6)

    # Cartesian RK4 two-body, converted back to elements, matches the analytic
    # Keplerian propagation (semi-major axis to sub-mm over the arc).
    coe_from_rv = np.asarray(pnt.cart_to_classical(rv_prop, GM))
    np.testing.assert_allclose(coe_from_rv[:, 0], coe_prop[:, 0], atol=1e-3)
    np.testing.assert_allclose(coe_from_rv[:, 1], coe_prop[:, 1], atol=1e-6)
    np.testing.assert_allclose(coe_from_rv[:, 2], coe_prop[:, 2], atol=1e-6)


def test_two_body_conserves_shape_elements():
    kep = pnt.KeplerianDynamics(GM)
    coe_prop = np.asarray(kep.propagate(COE0, TFS))
    # a, e, i, RAAN, argp are constant under unperturbed two-body motion.
    for col in range(5):
        np.testing.assert_allclose(coe_prop[:, col], COE0[col], atol=1e-6)


def test_cartesian_two_body_conserves_energy():
    rv0 = np.asarray(pnt.classical_to_cart(COE0, GM))
    cart = pnt.CartesianTwoBodyDynamics(GM)
    cart.set_time_step(5.0)
    rv_prop = np.asarray(cart.propagate(rv0, TFS))
    e0 = _specific_energy(rv0)
    for rv in rv_prop:
        assert abs(_specific_energy(rv) - e0) / abs(e0) < 1e-8


# ---------------------------------------------------------------------------
# J2 dynamics
# ---------------------------------------------------------------------------


def test_j2_zero_reduces_to_two_body():
    rv0 = np.asarray(pnt.classical_to_cart(COE0, GM))
    cart = pnt.CartesianTwoBodyDynamics(GM)
    cart.set_time_step(10.0)
    j2 = pnt.JToCartTwoBodyDynamics(GM, 0.0, R_BODY)
    j2.set_time_step(10.0)

    rv_two_body = np.asarray(cart.propagate(rv0, TFS))
    rv_j2_zero = np.asarray(j2.propagate(rv0, TFS))
    np.testing.assert_allclose(rv_j2_zero, rv_two_body, atol=1e-3)


def test_j2_nonzero_perturbs_orbit():
    rv0 = np.asarray(pnt.classical_to_cart(COE0, GM))
    cart = pnt.CartesianTwoBodyDynamics(GM)
    cart.set_time_step(10.0)
    j2 = pnt.JToCartTwoBodyDynamics(GM, 2.03e-4, R_BODY)  # lunar J2 ~2e-4
    j2.set_time_step(10.0)

    rv_two_body = np.asarray(cart.propagate(rv0, TFS))
    rv_j2 = np.asarray(j2.propagate(rv0, TFS))
    # A non-zero J2 measurably deflects the trajectory by the end of the arc.
    assert np.linalg.norm(rv_j2[-1, :3] - rv_two_body[-1, :3]) > 1.0
