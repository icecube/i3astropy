#!/usr/bin/env python

# SPDX-FileCopyrightText: © 2026 the i3astropy contributors
#
# SPDX-License-Identifier: BSD-2-Clause

"""Test Sidereal Functions in i3astropy."""

import contextlib
import sys

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import Longitude
from astropy.time import Time, TimeDelta

from i3astropy import gmast, gmest, gmst, gmt

with contextlib.suppress(ImportError):
    from icecube import astro  # pylint: disable=E0611
    from icecube.dataclasses import I3Time

approx = pytest.approx
DAYS_IN_YEAR = 365.2425
funcs = [
    (gmast, DAYS_IN_YEAR - 1, 1e-7),
    (gmt, DAYS_IN_YEAR, 1e-17),
    (gmst, DAYS_IN_YEAR + 1, 1e-7),
    (gmest, DAYS_IN_YEAR + 2, 1e-7),
]


def crossings(angles):
    """Calculate the number of times a function crosses 12."""
    anti = angles.hourangle - 12
    signchange = (anti[:-1] * anti[1:]) < 0
    positive = anti[:-1] < anti[1:]
    return np.logical_and(signchange, positive).sum()


@pytest.fixture
def setup_crossings():
    """Set up Crossings Test."""
    time0 = Time("2022-9-22 0:0:0", scale="utc")  # note: 4 years = 1461 days
    years = 4
    times = time0 + TimeDelta(
        np.linspace(0.03, years * 365.25 - 0.03, int(years * 365.25 * 10)),
        format="jd",
    )
    return years, times


@pytest.mark.parametrize(("func", "days", "_"), funcs, ids=[x[0].__name__ for x in funcs])
def test_crossings(setup_crossings, func, days, _):  # pylint: disable=W0621
    """Test the frequency of the function by testing the number of times the function crosses noon."""
    years, times = setup_crossings
    data = func(times)
    crossings_in_year = (crossings(data)) / years
    assert data[0] < 1 * u.hourangle
    assert data[-1] > 23 * u.hourangle
    assert crossings_in_year == days + 0.0075


@pytest.fixture
def setup_slope():
    """Set up the data for slope test."""
    time0 = Time("2022-9-20 12:0:0", scale="utc")
    times = time0 + TimeDelta(
        np.linspace(0, 2, 8192),
        format="jd",
    )
    return time0, times


@pytest.mark.parametrize(("func", "days", "tol"), funcs, ids=[x[0].__name__ for x in funcs])
def test_slope(setup_slope, func, days, tol):  # pylint: disable=W0621
    """Test the frequency of a function by measuring the slope in a single day."""
    time0, times = setup_slope
    rotation = func(times)
    indices = np.argwhere(np.diff(rotation) < 0)[:2, 0]
    indices[0] += 1
    rate = np.diff(rotation[indices].cycle)[0] / np.diff(times[indices].mjd)[0]

    assert isinstance(rotation, Longitude)
    assert rotation.unit == u.hourangle
    # make sure that the different timescales all almost match on september 21
    assert times[indices[0]].mjd == approx(int(time0.mjd) + 1, 1e-7)
    assert times[indices[-1]].mjd == approx(int(time0.mjd) + 2, 1e-7)
    # # make sure that the period matches what it is defined to be
    assert 365.2425 * rate == approx(days, tol)  # , delta=tol)


@pytest.fixture
def setup_fft():
    """Set up for samples for test with Fast Fourier Transform."""
    duration = 36524 // 2
    samples_per_day = 10
    samples = int(duration * samples_per_day)
    reftime = Time("1975-1-1 00:00:00", scale="utc")
    times = reftime + TimeDelta(np.linspace(0, duration, samples), format="jd")
    w = np.linspace(0.0, samples_per_day // 2, samples // 2) * DAYS_IN_YEAR
    return samples_per_day, times, w


@pytest.mark.parametrize(("func", "days", "_"), funcs, ids=[x[0].__name__ for x in funcs])
def test_fft(setup_fft, func, days, _):  # pylint: disable=W0621
    """Test Sidereal functions using fast fourier tranfrom to make sure the period is correct."""
    samples_per_day, times, w = setup_fft
    y = func(times).hour
    f = np.fft.rfft(y)[1:]
    i = np.abs(f).argmax()
    freq = w[i]
    assert freq == approx(days, samples_per_day / DAYS_IN_YEAR)


i3funcs = [
    gmast,
    gmst,
    gmest,
]


@pytest.mark.skipif("astro" not in globals(), reason="Not in an icetray invironment")
@pytest.mark.skipif(True, reason="Wait till #4269")  # noqa: FBT003
@pytest.mark.parametrize("func", [gmast, gmst, gmest])
def test_icetray(func):
    """Test sidereal functions against the ones in IceTray's Astro."""
    i3func = getattr(astro, "I3Get" + func.__name__.upper())
    time0 = Time("2022-9-20 12:0:0", scale="utc")
    times = time0 + TimeDelta(
        np.linspace(0, 365, 3650),
        format="jd",
    )
    for t in times:
        daq_year, daq_time = t.i3time
        i3t = I3Time(int(daq_year), int(daq_time))
        assert i3func(i3t) == approx(func(t).value, abs=1e-4)


if __name__ == "__main__":
    sys.exit(pytest.main(["-v", __file__, *sys.argv[1:]]))
