"""
pylupnt.plasma — ionospheric/plasmaspheric electron density and ray-tracing.

C++ implementation: cpp/lupnt/environment/plasma/ (namespace pecsim::)
Python bindings:    python/bindings/py_plasma.cc (InitPlasma → _pylupnt)
Data files:         $PECSIMPY_BASE_PATH/data/iri/  and  data/kp/
"""

from .kp_loader import convert_to_csv, update_kp, update_kp_table

# Physical / frequency constants re-exported with short names.
# The _PLASMA suffix on the _pylupnt side avoids collisions with other LuPNT
# symbols; here we expose them without the suffix for pecsim/plasma usage.
try:
    from pylupnt._pylupnt import (  # noqa: F401
        RE_PLASMA as RE,
        C_PLASMA as C,
        PI_PLASMA as PI,
        SECS_DAY_PLASMA as SECS_DAY,
        RAD2DEG_PLASMA as RAD2DEG,
        DEG2RAD_PLASMA as DEG2RAD,
        TECU_PLASMA as TECU,
        GM_EARTH_PLASMA as GM_EARTH,
        freq_L1,
        freq_L2,
        freq_L5,
    )

    _CONST_NAMES = [
        "RE",
        "C",
        "PI",
        "SECS_DAY",
        "RAD2DEG",
        "DEG2RAD",
        "TECU",
        "GM_EARTH",
        "freq_L1",
        "freq_L2",
        "freq_L5",
    ]
except ImportError:
    _CONST_NAMES = []

__all__ = [
    "update_kp",
    "update_kp_table",
    "convert_to_csv",
    # Frequency constants (same objects as the core _pylupnt symbols).
    "freq_L1",
    "freq_L2",
    "freq_L5",
    # NOTE: the generic physics constants RE / C / PI / SECS_DAY / RAD2DEG /
    # DEG2RAD / TECU / GM_EARTH are deliberately NOT listed in __all__. They are
    # the km-based pecsim values (C_PLASMA, GM_EARTH_PLASMA, ...). When they were
    # in __all__, the top-level package __init__'s `_try_import` merged them into
    # the `pylupnt` namespace via globals().update, shadowing the SI core
    # constants (e.g. making `pylupnt.C == 299792.458` km/s instead of the SI
    # `299792458` m/s, and `pylupnt.GM_EARTH == 398600.4418` km³/s²). Access the
    # plasma values explicitly instead: `pylupnt.plasma.C`, `pylupnt.plasma.RE`.
]

# Re-export the pybind11 symbols that live in the top-level _pylupnt module
# so users can also do `from pylupnt.plasma import trace_ray`, etc.
try:
    from pylupnt._pylupnt import (  # noqa: F401
        # core
        DateTime,
        get_plasma_base_path,
        set_plasma_base_path,
        # time utilities
        datetime_to_itime,
        itime_to_datetime,
        mjd_to_datetime,
        datetime_to_mjd,
        gregorian_to_mjd,
        tj2000_to_mjd,
        mjd_to_tj2000,
        long_to_lt,
        lt_to_long,
        sm_to_geo,
        geo_to_sm,
        # IRI / GCPM
        IRI2007Option,
        IRI2020Option,
        set_iri_model,
        get_iri_model,
        set_iri2007_option,
        set_iri2020_option,
        get_kp_index,
        gcpm_v24,
        gcpm_v24_fortran,
        # electron-density backend selection + NeQuick-G solar config
        set_iono_model,
        get_iono_model,
        NeQuickAzMode,
        NeQuickSolarConfig,
        # orbit utilities
        Satellite,
        wrap2pi,
        coe2cart,
        cart2coe,
        propagate_coe,
        compute_vis,
        compute_min_altitude,
        mean2true,
        true2mean,
        ecc2true,
        mean2ecc,
        solve_lt,
        setup_gnss_constellation,
        # TEC / ray-trace
        RayTraceConfig,
        PathProfile,
        trace_ray,
        compute_ne,
        get_iono_params,
        compute_B,
        refractive_index_neB,
    )

    __all__ += [
        # core
        "DateTime",
        "get_plasma_base_path",
        "set_plasma_base_path",
        # time utilities
        "datetime_to_itime",
        "itime_to_datetime",
        "mjd_to_datetime",
        "datetime_to_mjd",
        "gregorian_to_mjd",
        "tj2000_to_mjd",
        "mjd_to_tj2000",
        "long_to_lt",
        "lt_to_long",
        "sm_to_geo",
        "geo_to_sm",
        # IRI / GCPM
        "IRI2007Option",
        "IRI2020Option",
        "set_iri_model",
        "get_iri_model",
        "set_iri2007_option",
        "set_iri2020_option",
        "get_kp_index",
        "gcpm_v24",
        "gcpm_v24_fortran",
        # electron-density backend selection + NeQuick-G solar config
        "set_iono_model",
        "get_iono_model",
        "NeQuickAzMode",
        "NeQuickSolarConfig",
        # orbit utilities
        "Satellite",
        "wrap2pi",
        "coe2cart",
        "cart2coe",
        "propagate_coe",
        "compute_vis",
        "compute_min_altitude",
        "mean2true",
        "true2mean",
        "ecc2true",
        "mean2ecc",
        "solve_lt",
        "setup_gnss_constellation",
        # TEC / ray-trace
        "RayTraceConfig",
        "PathProfile",
        "trace_ray",
        "compute_ne",
        "get_iono_params",
        "compute_B",
        "refractive_index_neB",
    ]
    # Auto-configure the C++ base path when PECSIMPY_BASE_PATH is not set
    # (e.g. VSCode notebook kernels that run outside the pixi activation env).
    # Resolution order:
    #   1. PECSIMPY_BASE_PATH  (pixi activation — already handled by C++ init)
    #   2. LUPNT_DATA_PATH/plasma  (partial-env fallback)
    #   3. __file__-relative path  (works in any context; this file is always at
    #      <project_root>/python/pylupnt/plasma/__init__.py)
    import os as _os
    import pathlib as _pl

    if not get_plasma_base_path():
        _lupnt = _os.environ.get("LUPNT_DATA_PATH", "")
        _candidate = _os.path.join(_lupnt, "plasma") if _lupnt else ""
        if not (_candidate and _os.path.isdir(_candidate)):
            # Walk up: plasma/ → pylupnt/ → python/ → <project_root>
            _candidate = str(
                _pl.Path(__file__).resolve().parents[3] / "data" / "LuPNT_data" / "plasma"
            )
        if _os.path.isdir(_candidate):
            # Write into os.environ so that C++ get_base_path() — which
            # calls getenv("PECSIMPY_BASE_PATH") on *every* invocation —
            # always sees the correct value and never overrides the C++ cache.
            _os.environ["PECSIMPY_BASE_PATH"] = _candidate
            set_plasma_base_path(_candidate)

except ImportError:
    pass  # _pylupnt not yet built — that's OK at import time
