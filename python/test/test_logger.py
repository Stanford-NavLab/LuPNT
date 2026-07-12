"""Tests for pylupnt.core.logger (the Python logging helpers)."""

import pytest

from pylupnt.core import logger as lg


@pytest.mark.parametrize("seconds", [0.5, 65.0, 3661.0, 90061.0])
@pytest.mark.parametrize("short", [False, True])
def test_format_time_returns_str(seconds, short):
    out = lg.format_time(seconds, short=short)
    assert isinstance(out, str) and out


def test_logger_levels_do_not_raise():
    lg.Logger.set_log_level("DEBUG")
    lg.Logger.info("info message", name="test")
    lg.Logger.warning("warning message", name="test")
    lg.Logger.error("error message", name="test")
    lg.Logger.debug("debug message", name="test")
    lg.Logger.reset_start_time()


def test_logger_set_log_level_accepts_int_and_str():
    lg.Logger.set_log_level(20)
    lg.Logger.info("at INFO via int level", name="test")
    lg.Logger.set_log_level("WARNING")
    # below-threshold message must not raise even when filtered out
    lg.Logger.info("filtered out", name="test")
    lg.Logger.set_log_level("INFO")


def test_logger_log_with_time_sim():
    lg.Logger.set_log_level("INFO")
    lg.Logger.log("with sim time", name="test", time_sim=12.5)
