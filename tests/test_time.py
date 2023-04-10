#!/usr/bin/env python

# SPDX-FileCopyrightText: © 2022 the i3astropy contributors
#
# SPDX-License-Identifier: BSD-2-Clause

"""Tests related to I3Time in i3astropy."""

import contextlib

import i3astropy  # noqa: F401  pylint: disable=W0611
import numpy as np
import pytest
from astropy.time import Time
from astropy.time.core import ScaleValueError
from numpy.testing import assert_allclose, assert_equal

with contextlib.suppress(ImportError):
    from icecube.dataclasses import I3Time  # pylint: disable=E0611

DAQ_MS = int(1e7)
DAQ_SEC = int(1e10)
DAQ_MIN = 60 * DAQ_SEC
DAQ_HOUR = 60 * DAQ_MIN
DAQ_DAY = 24 * DAQ_HOUR

times = [
    (2015, 0, "2015-01-01 00:00:00"),
    (2015, DAQ_MS, "2015-01-01 00:00:00.001"),
    (2015, DAQ_SEC, "2015-01-01 00:00:01"),
    (2015, DAQ_MIN, "2015-01-01 00:01:00"),
    (2015, DAQ_HOUR, "2015-01-01 01:00:00"),
    (2015, DAQ_DAY, "2015-01-02 00:00:00"),
    (2015, 180 * DAQ_DAY, "2015-06-30 00:00:00"),
    (2015, 181 * DAQ_DAY - DAQ_SEC, "2015-06-30 23:59:59"),
    (2015, 181 * DAQ_DAY, "2015-06-30 23:59:60"),
    (2015, 181 * DAQ_DAY + DAQ_SEC, "2015-07-01 00:00:00"),
    (2015, 365 * DAQ_DAY, "2015-12-31 23:59:59"),
    (2016, 0, "2016-01-01 00:00:00"),
    (2016, DAQ_MS, "2016-01-01 00:00:00.001"),
    (2016, DAQ_SEC, "2016-01-01 00:00:01"),
    (2016, DAQ_MIN, "2016-01-01 00:01:00"),
    (2016, DAQ_HOUR, "2016-01-01 01:00:00"),
    (2016, DAQ_DAY, "2016-01-02 00:00:00"),
    (2016, 181 * DAQ_DAY, "2016-06-30 00:00:00"),
    (2016, 182 * DAQ_DAY, "2016-07-01 00:00:00"),
    (2016, 366 * DAQ_DAY - DAQ_SEC, "2016-12-31 23:59:59"),
    (2016, 366 * DAQ_DAY, "2016-12-31 23:59:60"),
    (2017, 0, "2017-01-01 00:00:00"),
    (2017, DAQ_MS, "2017-01-01 00:00:00.001"),
    (2017, DAQ_SEC, "2017-01-01 00:00:01"),
    (2017, DAQ_MIN, "2017-01-01 00:01:00"),
    (2017, DAQ_HOUR, "2017-01-01 01:00:00"),
    (2017, DAQ_DAY, "2017-01-02 00:00:00"),
    (2017, 180 * DAQ_DAY, "2017-06-30 00:00:00"),
    (2017, 181 * DAQ_DAY - DAQ_SEC, "2017-06-30 23:59:59"),
    (2017, 181 * DAQ_DAY, "2017-07-01 00:00:00"),
    (2017, 365 * DAQ_DAY - DAQ_SEC, "2017-12-31 23:59:59"),
]


def test_times(subtests):
    """Test I3Times match a list of iso times."""
    for daq_year, daq_time, iso in times:
        with subtests.test(daq_year=daq_year, daq_time=daq_time, iso=iso):
            i3t = Time(daq_year, daq_time, format="i3time")
            isot = Time(iso, format="iso", scale="utc")
            assert i3t.scale == "utc"
            assert np.all(i3t.isclose(isot))
            assert i3t.i3time.year == daq_year
            assert i3t.i3time.daq_time == daq_time
            assert isot.i3time.year == daq_year
            assert isot.i3time.daq_time == daq_time

    years, daqtime, iso = zip(*times)
    i3t = Time(years, daqtime, format="i3time")
    isot = Time(iso, format="iso")

    assert np.all(i3t.isclose(isot))
    assert_equal(i3t.i3time.year, years)
    assert_equal(i3t.i3time.daq_time, daqtime)
    assert_equal(isot.i3time.year, years)
    assert_equal(isot.i3time.daq_time, daqtime)

    with pytest.raises(ScaleValueError):
        Time(2020, 0, format="i3time", scale="tt")
    with pytest.raises(ValueError, match="2020.5"):
        Time(2020.5, 0, format="i3time")


@pytest.mark.skipif("I3Time" not in globals(), reason="Not in an icetray invironment")
def test_icetray(subtests):
    """Test i3astropy.I3Time matches dataclasses.I3Time."""
    for daq_year, daq_time, iso in times:
        with subtests.test(daq_year=daq_year, daq_time=daq_time, iso=iso):
            di3t = I3Time(daq_year, daq_time)
            ai3t = Time(daq_year, daq_time, format="i3time")
            isot = Time(iso, format="iso")
            assert str(di3t)[: len(iso)] == iso
            assert_allclose(di3t.mod_julian_day_double, ai3t.mjd, rtol=1e-9)
            assert_allclose(di3t.mod_julian_day_double, isot.mjd, rtol=1e-9)


if __name__ == "__main__":
    pytest.main()
