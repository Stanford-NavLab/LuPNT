import pylupnt as pnt
import numpy as np
import pytest
from ..gmat.utils import data
from test.gmat.utils import data

try:
    from .utils import gmat_helpers
    from .utils.gmat import gmat
except ImportError:
    from utils import gmat_helpers
    from utils.gmat import gmat


def convert_time_gmat(epoch, time_sys_from, time_sys_to):
    time_converter = gmat.TimeSystemConverter.Instance()
    return (
        time_converter.Convert(epoch, time_sys_from, time_sys_to) * gmat.SECS_PER_DAY
        - gmat.MJD_OF_J2000 * gmat.SECS_PER_DAY
    )


class TestSpiceInterface:
    def test_time_conversions(self):
        string = "2020/07/20 12:00:00.000 UTC"
        string_gmat = "20 Jul 2020 12:00:00.000"

        time_converter = gmat.TimeSystemConverter.Instance()

        # String UTC to TAI
        epoch = pnt.SpiceInterface.string_to_tai(string)
        epoch_gmat = gmat_helpers.convert_time(
            time_converter.ConvertGregorianToMjd(string_gmat),
            gmat.TimeSystemConverter.UTC,
            gmat.TimeSystemConverter.TAI,
        )

        assert np.isclose(epoch, epoch_gmat)

        # String UTC to TDB
        epoch = pnt.SpiceInterface.string_to_tdb(string)
        epoch_gmat = gmat_helpers.convert_time(
            time_converter.ConvertGregorianToMjd(string_gmat),
            gmat.TimeSystemConverter.UTC,
            gmat.TimeSystemConverter.TDB,
        )
        assert np.isclose(epoch, epoch_gmat)


if __name__ == "__main__":
    pytest.main([__file__])
