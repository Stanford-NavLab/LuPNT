"""Unit tests for the `pylupnt` Earth-orientation bindings (`get_eop_coverage`,
`get_eop_data`, `get_eop_file_data`, `load_eop_file_data`,
`load_eop_finals_file_data`).

These exercise the two EOP sources through their Python bindings: the bundled
IERS C04 final series and the IERS finals ("Bulletin A") reader, which is the
only one carrying predictions. Network-dependent loaders
(`load_latest_eop_*_from_iers`) are deliberately not exercised here.

Every test restores the bundled C04 table on the way out, since the table is
global process state and the rest of the suite is calibrated against it.
"""

import os

import numpy as np
import pylupnt as pnt
import pytest

C04_FILENAME = "eopc04_08.62-now"

# The repo's orekit data snapshot ships an IAU-1980 finals.all, used here so the
# finals reader can be tested without network access.
FINALS_PATH = os.path.join(
    pnt.get_data_path(),
    "..",
    "orekit-data-main",
    "Earth-Orientation-Parameters",
    "IAU-1980",
    "finals.all",
)

RAD_ARCSEC = pnt.RAD_ARCSEC


@pytest.fixture(autouse=True)
def restore_bundled_eop():
    """Reload the bundled C04 table after each test."""
    yield
    pnt.load_eop_file_data(pnt.get_file_path(C04_FILENAME), True)


def _load_c04():
    pnt.load_eop_file_data(pnt.get_file_path(C04_FILENAME), True)


def _load_finals():
    if not os.path.exists(FINALS_PATH):
        pytest.skip("orekit finals.all snapshot not present")
    pnt.load_eop_finals_file_data(FINALS_PATH, True)


def test_c04_coverage_and_source():
    _load_c04()
    cov = pnt.get_eop_coverage()
    assert cov["source"] == pnt.EopSource.C04
    assert cov["mjd_first"] == pytest.approx(37665.0)  # 1962-01-01
    assert cov["mjd_last"] > cov["mjd_first"]
    # C04 marks non-final rows with a 0.999 formal-error sentinel, so the
    # measured span ends at or before the table end.
    assert cov["mjd_last_measured"] <= cov["mjd_last"]


def test_eop_data_keys_and_units():
    _load_c04()
    cov = pnt.get_eop_coverage()
    data = pnt.get_eop_data(cov["mjd_last_measured"] - 100.0)
    assert set(data) == {
        "x_pole",
        "y_pole",
        "ut1_utc",
        "lod",
        "dpsi",
        "deps",
        "sigma_x_pole",
        "sigma_y_pole",
        "sigma_ut1_utc",
        "dX",
        "dY",
    }
    # C04 carries dpsi/deps, not CIP offsets, so dX/dY are zero for this source.
    assert data["dX"] == 0.0 and data["dY"] == 0.0
    # Polar motion is well under an arcsecond; UT1-UTC is kept under a second by
    # leap seconds; LOD is sub-millisecond.
    assert abs(data["x_pole"] / RAD_ARCSEC) < 1.0
    assert abs(data["y_pole"] / RAD_ARCSEC) < 1.0
    assert abs(data["ut1_utc"]) < 1.0
    assert abs(data["lod"]) < 1e-2
    # The sigmas are formal errors, and on final rows they are far smaller than
    # the values themselves.
    assert 0.0 < data["sigma_x_pole"] / RAD_ARCSEC < 1e-2
    assert 0.0 < data["sigma_ut1_utc"] < 1e-3


def test_out_of_range_epochs_are_clamped_not_extrapolated():
    _load_c04()
    cov = pnt.get_eop_coverage()
    at_end = pnt.get_eop_data(cov["mjd_last"])
    beyond = pnt.get_eop_data(cov["mjd_last"] + 500.0)
    way_beyond = pnt.get_eop_data(cov["mjd_last"] + 5000.0)
    # Held constant: both out-of-range requests return the same endpoint values.
    for key in ("x_pole", "y_pole", "ut1_utc", "lod"):
        assert beyond[key] == pytest.approx(at_end[key], abs=1e-12)
        assert way_beyond[key] == pytest.approx(at_end[key], abs=1e-12)


def test_finals_reader_reaches_further_than_c04():
    _load_c04()
    c04_last = pnt.get_eop_coverage()["mjd_last"]

    _load_finals()
    cov = pnt.get_eop_coverage()
    assert cov["source"] == pnt.EopSource.Finals
    assert cov["mjd_first"] == pytest.approx(41684.0)  # 1973-01-02
    # The point of the finals product: it reaches past the C04 cutoff and keeps
    # going into a predicted span.
    assert cov["mjd_last_measured"] > c04_last
    assert cov["mjd_last"] > cov["mjd_last_measured"]


def test_finals_prediction_flags_are_a_trailing_block():
    _load_finals()
    fd = pnt.get_eop_file_data()
    cov = pnt.get_eop_coverage()
    mjds, is_pred = fd["mjds_utc"], fd["is_prediction"].astype(bool)

    assert fd["source"] == pnt.EopSource.Finals
    assert len(is_pred) == len(mjds)
    assert np.all(np.diff(mjds) > 0)
    assert is_pred.any() and not is_pred.all()
    # Predictions are a single trailing block, delimited by mjd_last_measured.
    assert np.all(np.diff(is_pred.astype(int)) >= 0)
    np.testing.assert_array_equal(is_pred, mjds > cov["mjd_last_measured"])


def test_finals_values_are_physical_including_predictions():
    _load_finals()
    fd = pnt.get_eop_file_data()
    assert np.all(np.abs(fd["x"]) < 1.0)  # ["]
    assert np.all(np.abs(fd["y"]) < 1.0)
    assert np.all(np.abs(fd["ut1_utc"]) < 1.0)  # [s]
    assert np.all(np.abs(fd["lod"]) < 1e-2)


def test_finals_agrees_with_c04_where_both_are_final():
    """MJD 59000 (2020-05-31) is well inside both series' final spans."""
    mjd = 59000.0
    _load_finals()
    finals = pnt.get_eop_data(mjd)
    _load_c04()
    c04 = pnt.get_eop_data(mjd)

    # Tolerances are the genuine C04-vs-Bulletin-A solution difference.
    assert (finals["x_pole"] - c04["x_pole"]) / RAD_ARCSEC == pytest.approx(0.0, abs=5e-3)
    assert (finals["y_pole"] - c04["y_pole"]) / RAD_ARCSEC == pytest.approx(0.0, abs=5e-3)
    assert finals["ut1_utc"] - c04["ut1_utc"] == pytest.approx(0.0, abs=1e-4)
    assert finals["lod"] - c04["lod"] == pytest.approx(0.0, abs=1e-4)


def test_set_eop_source_defaults_to_c04():
    _load_c04()
    assert pnt.get_eop_source() == pnt.EopSource.C04
    assert pnt.get_eop_coverage()["source"] == pnt.EopSource.C04


def test_set_eop_source_switches_an_already_loaded_table():
    """Order independence: selecting a source with a table already loaded must
    switch immediately, not defer (which would amount to ignoring the call)."""
    _load_c04()
    if not os.path.exists(FINALS_PATH):
        pytest.skip("orekit finals.all snapshot not present")

    pnt.set_eop_source(pnt.EopSource.Finals, FINALS_PATH)
    assert pnt.get_eop_source() == pnt.EopSource.Finals
    assert pnt.get_eop_coverage()["source"] == pnt.EopSource.Finals

    pnt.set_eop_source(pnt.EopSource.C04)
    assert pnt.get_eop_coverage()["source"] == pnt.EopSource.C04


def test_set_eop_source_rejects_a_missing_file_and_keeps_the_previous_choice():
    _load_c04()
    with pytest.raises(RuntimeError):
        pnt.set_eop_source(pnt.EopSource.Finals, "/nonexistent/finals.all")
    assert pnt.get_eop_source() == pnt.EopSource.C04
    assert pnt.get_eop_coverage()["source"] == pnt.EopSource.C04


def test_non_forcing_load_is_dropped_when_a_table_is_present():
    """The trap set_eop_source routes around: load_*(force=False) after anything
    has touched EOP is a no-op, leaving the previous product in place."""
    _load_c04()
    if not os.path.exists(FINALS_PATH):
        pytest.skip("orekit finals.all snapshot not present")

    pnt.load_eop_finals_file_data(FINALS_PATH)  # force defaults to False
    assert pnt.get_eop_coverage()["source"] == pnt.EopSource.C04

    pnt.load_eop_finals_file_data(FINALS_PATH, True)
    assert pnt.get_eop_coverage()["source"] == pnt.EopSource.Finals


def test_finals_predicted_sigmas_grow_with_lead_time():
    """Bulletin A publishes its own predicted formal errors, which grow with
    lead time -- the first-cut EOP prediction error budget."""
    _load_finals()
    cov = pnt.get_eop_coverage()
    horizon = cov["mjd_last"] - cov["mjd_last_measured"]
    leads = [lead for lead in (1.0, 30.0, 90.0, 180.0) if lead < horizon]
    sigmas = [pnt.get_eop_data(cov["mjd_last_measured"] + lead)["sigma_ut1_utc"] for lead in leads]

    assert np.all(np.diff(sigmas) > 0), f"sigma_ut1_utc not increasing: {sigmas}"
    # Sanity anchors on the IERS-published magnitudes: sub-ms at a day, several
    # ms at a month.
    assert sigmas[0] < 1e-3
    assert 1e-3 < sigmas[leads.index(30.0)] < 2e-2


# ---------------------------------------------------------------------------
# Step 1 / Step 2a: celestial-pole offsets and the perturbation injector.
# ---------------------------------------------------------------------------

FINALS2000A_PATH = os.path.join(
    pnt.get_data_path(),
    "..",
    "orekit-data-main",
    "Earth-Orientation-Parameters",
    "IAU-2000",
    "finals2000A.all",
)


def test_iau2000_source_populates_cip_offsets():
    """finals2000A carries dX/dY where finals.all carries dpsi/deps; the reader
    must file them under whichever the caller declared."""
    if not os.path.exists(FINALS2000A_PATH):
        pytest.skip("orekit finals2000A.all snapshot not present")
    pnt.load_eop_finals_file_data(FINALS2000A_PATH, True, pnt.EopNutation.Iau2000)
    fd = pnt.get_eop_file_data()
    assert fd["nutation"] == pnt.EopNutation.Iau2000
    assert np.any(fd["dX"] != 0.0)
    assert np.all(fd["dpsi"] == 0.0)  # the unused pair is zero-filled, not left stale
    # CIP offsets are ~0.3 mas rms; anything much larger means dpsi/deps got mis-filed.
    assert np.percentile(np.abs(fd["dX"]), 99) < 5.0  # ["] -- generous, catches a 1000x mix-up
    assert pnt.eop_has_celestial_pole_offsets()

    data = pnt.get_eop_data(59000.0)
    assert abs(data["dX"] / RAD_ARCSEC * 1e3) < 5.0  # mas
    assert abs(data["dY"] / RAD_ARCSEC * 1e3) < 5.0


def test_perturbation_is_additive_and_clearable():
    _load_c04()
    pnt.clear_eop_perturbation()
    assert not pnt.eop_has_celestial_pole_offsets()
    nominal = pnt.get_eop_data(59000.0)

    p = pnt.EopPerturbation()
    p.dx_pole, p.dy_pole, p.dut1, p.dlod = 1e-8, -2e-8, 3e-3, 4e-6
    pnt.set_eop_perturbation(p)
    got = pnt.get_eop_data(59000.0)
    assert got["x_pole"] - nominal["x_pole"] == pytest.approx(1e-8, abs=1e-14)
    assert got["y_pole"] - nominal["y_pole"] == pytest.approx(-2e-8, abs=1e-14)
    assert got["ut1_utc"] - nominal["ut1_utc"] == pytest.approx(3e-3, abs=1e-9)
    assert got["lod"] - nominal["lod"] == pytest.approx(4e-6, abs=1e-12)

    pnt.clear_eop_perturbation()
    restored = pnt.get_eop_data(59000.0)
    assert restored["x_pole"] == pytest.approx(nominal["x_pole"], abs=1e-15)
    assert restored["ut1_utc"] == pytest.approx(nominal["ut1_utc"], abs=1e-12)


def test_perturbation_reaches_the_frame_rotation_with_the_expected_sign():
    """The small-angle convention, pinned by central difference rather than asserted:
    delta r_gcrf = R_itrf->gcrf @ (eps x r_itrf), eps = (-dy_pole, -dx_pole, OMEGA*dut1)."""
    _load_c04()
    omega = 2 * np.pi * 1.00273781191135448 / 86400.0
    t_tdb = pnt.convert_time(pnt.mjd_to_time(59000.0), pnt.Time.UTC, pnt.Time.TDB)
    r_itrf = np.array([4197160.8, 815845.4, 4716876.3])

    def rot(**kw):
        pnt.clear_eop_perturbation()
        if kw:
            p = pnt.EopPerturbation()
            for k, v in kw.items():
                setattr(p, k, v)
            pnt.set_eop_perturbation(p)
        R, _ = pnt.get_frame_rotation_translation(t_tdb, pnt.Frame.ITRF, pnt.Frame.GCRF)
        pnt.clear_eop_perturbation()
        return np.asarray(R, float)

    R = rot()
    # Steps must clear double-precision noise: epochs are ~1.9e9 s since J2000, so a
    # 1e-6 s UT1 step is pure roundoff. See step2a_sensitivity.py.
    for name, step, axis in [
        ("dx_pole", 1e-8, np.array([0.0, -1.0, 0.0])),
        ("dy_pole", 1e-8, np.array([-1.0, 0.0, 0.0])),
        ("dut1", 1e-2, np.array([0.0, 0.0, omega])),
    ]:
        num = (rot(**{name: step}) @ r_itrf - rot(**{name: -step}) @ r_itrf) / (2 * step)
        ana = R @ np.cross(axis, r_itrf)
        rel = np.linalg.norm(num - ana) / np.linalg.norm(ana)
        assert rel < 1e-4, f"{name}: relative error {rel:.2e}"

    pnt.clear_eop_perturbation()
