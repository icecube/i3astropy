#!/usr/bin/env python

# SPDX-FileCopyrightText: © 2022 the i3astropy contributors
#
# SPDX-License-Identifier: BSD-2-Clause

"""Test Coordinate transforms in i3astropy."""

import contextlib
import sys

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import ICRS, Angle, SkyCoord, get_body, get_sun
from astropy.time import Time
from astropy.units import day, deg, hour

from i3astropy import I3Dir

with contextlib.suppress(ImportError):
    from icecube import astro
approx = pytest.approx


def test_j2000_to_i3dir():
    """Test conversions from J2000 to I3Direction."""
    obs_time = Time("2020-03-22 12:00:00", format="iso", scale="utc")

    # The first point of Ares should be grid north at noon on the vernal equinox
    i3fpa = SkyCoord(ra=0 * u.deg, dec=0 * u.deg, frame="icrs", obstime=obs_time).transform_to(I3Dir())
    assert i3fpa.zen.degree == approx(90, abs=0.15)
    assert i3fpa.az.degree == approx(90, abs=0.2)

    # RA = 90 should be Grid East
    i3ra90 = SkyCoord(ra=90 * u.deg, dec=0 * u.deg, frame="icrs", obstime=obs_time).transform_to(
        I3Dir(),
    )
    assert i3ra90.zen.degree == approx(90, abs=0.003)
    assert i3ra90.az.degree == approx(0, abs=0.2)

    # RA =180 should be grid south
    i3ra180 = SkyCoord(ra=180 * u.deg, dec=0 * u.deg, frame="icrs", obstime=obs_time).transform_to(
        I3Dir(),
    )
    assert i3ra180.zen.degree == approx(90, abs=0.15)
    assert i3ra180.az.degree == approx(270, abs=0.2)

    # RA =270 should be grid west
    i3ra270 = SkyCoord(ra=270 * u.deg, dec=0 * u.deg, frame="icrs", obstime=obs_time).transform_to(
        I3Dir(),
    )
    assert i3ra270.zen.degree == approx(90, abs=0.01)
    assert i3ra270.az.degree == approx(180, abs=0.2)

    # celestial north pole should be nadir
    i3np = SkyCoord(ra=0 * u.deg, dec=+90 * u.deg, frame="icrs", obstime=obs_time).transform_to(I3Dir())
    assert i3np.zen.degree == approx(180, abs=0.15)

    # celestial south pole should be zenith
    i3sp = SkyCoord(ra=0 * u.deg, dec=-90 * u.deg, frame="icrs", obstime=obs_time).transform_to(I3Dir())
    assert i3sp.zen.degree == approx(0, abs=0.15)


def test_j2000_to_i3dir_array():
    """Test conversions from J2000 to I3Direction with arrays."""
    obs_time = Time("2020-03-22 12:00:00", format="iso", scale="utc")

    ras = [0, 90, 180, 270, 0, 0] * u.deg
    dec = [+0, +0, +0, +0, +90, -90] * u.deg
    azi = [90, 0, 270, 180]
    zen = [90, 90, 90, 90, 180, 0]

    i3dir = SkyCoord(ra=ras, dec=dec, frame="icrs", obstime=obs_time).transform_to(I3Dir())
    assert i3dir.az.degree[:4] == approx(azi, abs=0.2)
    assert i3dir.zen.degree == approx(zen, abs=0.2)

    i3dir = ICRS(ra=ras, dec=dec).transform_to(I3Dir(obstime=obs_time))
    assert i3dir.az.degree[:4] == approx(azi, abs=0.2)
    assert i3dir.zen.degree == approx(zen, abs=0.2)

    obs_time1 = obs_time + range(25) * hour
    i3dir = SkyCoord(ra=0 * deg, dec=0 * deg, obstime=obs_time1).transform_to(I3Dir())
    azimuth = Angle(np.arange(90, 452, 15 + 1 / 24), unit=deg)
    azimuth = azimuth.wrap_at(360 * deg).degree
    assert i3dir.az.degree == approx(azimuth, abs=0.2)
    assert i3dir.zen.degree == approx(90, abs=0.2)

    obs_time2 = obs_time + range(366) * day
    i3dir = SkyCoord(ra=0 * deg, dec=0 * deg).transform_to(I3Dir(obstime=obs_time2))
    azimuth = Angle(np.linspace(90, 450, 366), unit=deg)
    azimuth = azimuth.wrap_at(360 * deg).degree
    assert i3dir.az.degree == approx(azimuth, abs=0.2)
    assert i3dir.zen.degree == approx(90, abs=0.2)


def test_i3dir_to_j2000():
    """Conversions from I3Direction to J2000."""
    obs_time = Time("2020-03-22 12:00:00", format="iso", scale="utc")

    # grid east should be ra=90
    grid_east = I3Dir(zen=90 * u.deg, az=0 * u.deg, obstime=obs_time).transform_to(ICRS())
    assert grid_east.ra.degree == approx(90, abs=0.2)
    assert grid_east.dec.degree == approx(0, abs=0.01)

    # grid north should be Zero Point of Ares
    grid_north = I3Dir(zen=90 * u.deg, az=90 * u.deg, obstime=obs_time).transform_to(ICRS())
    assert grid_north.ra.degree == approx(0, abs=0.2)
    assert grid_north.dec.degree == approx(0, abs=0.15)

    # grid west should be ra=270
    grid_west = I3Dir(zen=90 * u.deg, az=180 * u.deg, obstime=obs_time).transform_to(ICRS())
    assert grid_west.ra.degree == approx(270, abs=0.2)
    assert grid_west.dec.degree == approx(0, abs=0.02)

    # grid south should be ra=180
    grid_south = I3Dir(zen=90 * u.deg, az=270 * u.deg, obstime=obs_time).transform_to(ICRS())
    assert grid_south.ra.degree == approx(180, abs=0.2)
    assert grid_south.dec.degree == approx(0, abs=0.15)

    # zenith should be celestial south pole
    zenith = I3Dir(zen=0 * u.deg, az=0 * u.deg, obstime=obs_time).transform_to(ICRS())
    assert zenith.dec.degree == approx(-90, abs=0.15)

    # nadir should be celestial south pole
    nadir = I3Dir(zen=180 * u.deg, az=0 * u.deg, obstime=obs_time).transform_to(ICRS())
    assert nadir.dec.degree == approx(+90, abs=0.15)


def test_sun():
    """Conversions from the sun to I3Direction."""
    # times when the Equation of time is stationary
    assert get_sun(Time("2020-04-15 12:00")).transform_to(I3Dir()).az.degree == approx(90, abs=0.02)
    assert get_sun(Time("2020-06-13 12:00")).transform_to(I3Dir()).az.degree == approx(90, abs=0.05)
    assert get_sun(Time("2020-09-01 12:00")).transform_to(I3Dir()).az.degree == approx(90, abs=0.04)
    assert get_sun(Time("2020-12-24 12:00")).transform_to(I3Dir()).az.degree == approx(90, abs=0.05)

    # times when the Equation of time is maximum/minimum
    assert get_sun(Time("2020-02-11 12:14:15")).transform_to(I3Dir()).az.degree == approx(90, abs=0.02)
    assert get_sun(Time("2020-05-14 11:56:19")).transform_to(I3Dir()).az.degree == approx(90, abs=0.02)
    assert get_sun(Time("2020-07-26 12:06:36")).transform_to(I3Dir()).az.degree == approx(90, abs=0.02)
    assert get_sun(Time("2020-11-03 11:43:35")).transform_to(I3Dir()).az.degree == approx(90, abs=0.02)

    assert get_sun(Time("2020-03-20 03:50")).transform_to(I3Dir()).zen.degree == approx(90, abs=0.02)
    assert get_sun(Time("2020-06-20 21:43")).transform_to(I3Dir()).zen.degree == approx(113.44, abs=1)
    assert get_sun(Time("2020-09-22 13:31")).transform_to(I3Dir()).zen.degree == approx(90, abs=0.02)
    assert get_sun(Time("2020-12-21 10:03")).transform_to(I3Dir()).zen.degree == approx(66.56, abs=0.02)


def test_sun_array():
    """Conversions from the sun to I3Direction with arrays."""
    ref_time = Time("2020-01-01 12:00")
    day_offsets = np.arange(366)
    obs_time = ref_time + day_offsets * day
    sun1 = get_sun(obs_time).transform_to(I3Dir())
    d = 6.240_040_77 + 0.017_201_97 * (365.25 * (ref_time.ymdhms.year - 2000) + day_offsets)
    ref1 = I3Dir(
        zen=(90 + 23.44 * np.sin((day_offsets - 80) / len(day_offsets) * 2 * np.pi)) * deg,
        az=(90 + (-7.659 * np.sin(d) + 9.863 * np.sin(2 * d + 3.5932)) / 4) * u.deg,
    )

    assert sun1.zen.degree == approx(ref1.zen.degree, rel=0.02)
    assert sun1.az.degree == approx(ref1.az.degree, rel=0.02)
    assert sun1.separation(ref1).degree == approx(0, abs=1)

    ref_time = Time("2020-03-20 00:00")
    day_inc = np.linspace(0, 1, 1441)
    obs_time2 = ref_time + day_inc * day
    sun2 = get_sun(obs_time2).transform_to(I3Dir())
    ref2 = I3Dir(zen=90 * deg, az=(268.2 + day_inc * 360) * deg)
    assert sun2.zen.degree == approx(ref2.zen.degree, rel=0.004)
    assert (sun2.az - ref2.az).wrap_at(180 * deg).degree == approx(0, abs=0.1)
    assert sun2.separation(ref2).degree == approx(0, abs=0.4)


@pytest.mark.skipif("astro" not in globals(), reason="Not in an icetray invironment")
def test_icetray():
    # pylint: disable=R0914
    """Comparisons with IceTray's astro module."""
    obs_time = Time("2020-03-20 0:00") + np.arange(366) * u.day
    rng = np.random.default_rng(seed=0)
    zen = 1 * np.pi * rng.random(len(obs_time))
    azi = 2 * np.pi * rng.random(len(obs_time))

    i3dir = I3Dir(zen=zen * u.radian, az=azi * u.radian, obstime=obs_time)
    equa = i3dir.transform_to(ICRS())
    ras, dec = astro.dir_to_equa(zen, azi, obs_time.mjd)

    assert equa.ra.radian == approx(ras, abs=1e-3)
    assert equa.dec.radian == approx(dec, abs=1e-5)

    crab = SkyCoord.from_name("Crab")
    i3crab = crab.transform_to(I3Dir(obstime=obs_time))
    zenith, azimuth = astro.equa_to_dir(crab.ra.radian, crab.dec.radian, obs_time.mjd)
    assert zenith == approx(i3crab.zen.radian, abs=1e-5)
    assert azimuth == approx(i3crab.az.radian, abs=2e-5)

    i3sun = get_sun(obs_time).transform_to(I3Dir())
    sun_zen, sun_azi = astro.sun_dir(obs_time.mjd)
    assert sun_zen == approx(i3sun.zen.radian, abs=2e-5)
    assert sun_azi == approx(i3sun.az.radian, abs=1e-4)

    i3moon = get_body("moon", obs_time).transform_to(I3Dir())
    moon_zen, moon_azi = astro.moon_dir(obs_time.mjd)
    assert moon_zen == approx(i3moon.zen.radian, abs=1e-4)
    assert moon_azi == approx(i3moon.az.radian, abs=1e-4)


if __name__ == "__main__":
    sys.exit(pytest.main(["-v", __file__, *sys.argv[1:]]))
