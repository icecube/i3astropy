#!/usr/bin/env python

# Copyright (c) 2022, The IceCube Collaboration, All Rights Reserved
# SPDX-License-Identifier:  GPL-3.0-or-later
# Author: Kevin Meagher

import unittest

import numpy as np
from astropy import units as u
from astropy.coordinates import ICRS, Angle, SkyCoord, get_moon, get_sun
from astropy.time import Time
from astropy.units import day, deg, hour
from numpy.testing import assert_allclose

from i3astropy import I3Dir

try:
    from icecube import astro
except ImportError:
    pass


class TestI3AstroPy(unittest.TestCase):
    def test_j2000_to_i3dir(self):

        t = Time("2020-03-22 12:00:00", format="iso", scale="utc")

        # The first point of Ares should be grid north at noon on the vernal equinox
        i3fpa = SkyCoord(ra=0 * u.deg, dec=0 * u.deg, frame="icrs", obstime=t).transform_to(I3Dir())
        self.assertAlmostEqual(i3fpa.zen.degree, 90, 0)
        self.assertAlmostEqual(i3fpa.az.degree, 90, 0)

        # RA = 90 should be Grid East
        i3ra90 = SkyCoord(ra=90 * u.deg, dec=0 * u.deg, frame="icrs", obstime=t).transform_to(I3Dir())
        self.assertAlmostEqual(i3ra90.zen.degree, 90, 0)
        self.assertAlmostEqual(i3ra90.az.degree, 0, 0)

        # RA =180 should be grid south
        i3ra180 = SkyCoord(ra=180 * u.deg, dec=0 * u.deg, frame="icrs", obstime=t).transform_to(I3Dir())
        self.assertAlmostEqual(i3ra180.zen.degree, 90, 0)
        self.assertAlmostEqual(i3ra180.az.degree, 270, 0)

        # RA =270 should be grid west
        i3ra270 = SkyCoord(ra=270 * u.deg, dec=0 * u.deg, frame="icrs", obstime=t).transform_to(I3Dir())
        self.assertAlmostEqual(i3ra270.zen.degree, 90, 0)
        self.assertAlmostEqual(i3ra270.az.degree, 180, 0)

        # celestial north pole should be nadir
        i3np = SkyCoord(ra=0 * u.deg, dec=+90 * u.deg, frame="icrs", obstime=t).transform_to(I3Dir())
        self.assertAlmostEqual(i3np.zen.degree, 180, 0)

        # celestial south pole should be zenith
        i3sp = SkyCoord(ra=0 * u.deg, dec=-90 * u.deg, frame="icrs", obstime=t).transform_to(I3Dir())
        self.assertAlmostEqual(i3sp.zen.degree, 0, 0)

    def test_j2000_to_i3dir_array(self):

        t = Time("2020-03-22 12:00:00", format="iso", scale="utc")

        ra = [0, 90, 180, 270, 0, 0] * u.deg
        dec = [+0, +0, +0, +0, +90, -90] * u.deg
        azi = [90, 0, 270, 180]
        zen = [90, 90, 90, 90, 180, 0]

        i3dir = SkyCoord(ra=ra, dec=dec, frame="icrs", obstime=t).transform_to(I3Dir())
        assert_allclose(i3dir.az.degree[:4], azi, atol=0.2)
        assert_allclose(i3dir.zen.degree, zen, atol=0.2)

        i3dir = ICRS(ra=ra, dec=dec).transform_to(I3Dir(obstime=t))
        assert_allclose(i3dir.az.degree[:4], azi, atol=0.2)
        assert_allclose(i3dir.zen.degree, zen, atol=0.2)

        t1 = t + range(25) * hour
        i3dir = SkyCoord(ra=0 * deg, dec=0 * deg, obstime=t1).transform_to(I3Dir())
        x = Angle(np.arange(90, 452, 15 + 1 / 24), unit=deg)
        x = x.wrap_at(360 * deg).degree
        assert_allclose(i3dir.az.degree, x, atol=0.2)
        assert_allclose(i3dir.zen.degree, 90, atol=0.2)

        t2 = t + range(366) * day
        i3dir = SkyCoord(ra=0 * deg, dec=0 * deg).transform_to(I3Dir(obstime=t2))
        x = Angle(np.linspace(90, 450, 366), unit=deg)
        x = x.wrap_at(360 * deg).degree
        assert_allclose(i3dir.az.degree, x, atol=0.2)
        assert_allclose(i3dir.zen.degree, 90, atol=0.2)

    def test_i3dir_to_j2000(self):

        t = Time("2020-03-22 12:00:00", format="iso", scale="utc")

        # grid east should be ra=90
        ge = I3Dir(zen=90 * u.deg, az=0 * u.deg, obstime=t).transform_to(ICRS())
        self.assertAlmostEqual(ge.ra.degree, 90, 0)
        self.assertAlmostEqual(ge.dec.degree, 0, 0)

        # grid north should be Zero Point of Ares
        ge = I3Dir(zen=90 * u.deg, az=90 * u.deg, obstime=t).transform_to(ICRS())
        self.assertAlmostEqual(ge.ra.degree, 0, 0)
        self.assertAlmostEqual(ge.dec.degree, 0, 0)

        # grid west should be ra=270
        ge = I3Dir(zen=90 * u.deg, az=180 * u.deg, obstime=t).transform_to(ICRS())
        self.assertAlmostEqual(ge.ra.degree, 270, 0)
        self.assertAlmostEqual(ge.dec.degree, 0, 0)

        # grid south should be ra=180
        ge = I3Dir(zen=90 * u.deg, az=270 * u.deg, obstime=t).transform_to(ICRS())
        self.assertAlmostEqual(ge.ra.degree, 180, 0)
        self.assertAlmostEqual(ge.dec.degree, 0, 0)

        # zenith should be celestial south pole
        ge = I3Dir(zen=0 * u.deg, az=0 * u.deg, obstime=t).transform_to(ICRS())
        self.assertAlmostEqual(ge.dec.degree, -90, 0)

        # nadir should be celestial south pole
        ge = I3Dir(zen=180 * u.deg, az=0 * u.deg, obstime=t).transform_to(ICRS())
        self.assertAlmostEqual(ge.dec.degree, +90, 0)

    def test_sun(self):
        # times when the Equation of time is stationary
        self.assertAlmostEqual(get_sun(Time("2020-04-15 12:00")).transform_to(I3Dir()).az.degree, 90, 1)
        self.assertAlmostEqual(get_sun(Time("2020-06-13 12:00")).transform_to(I3Dir()).az.degree, 90, 1)
        self.assertAlmostEqual(get_sun(Time("2020-09-01 12:00")).transform_to(I3Dir()).az.degree, 90, 1)
        self.assertAlmostEqual(get_sun(Time("2020-12-24 12:00")).transform_to(I3Dir()).az.degree, 90, 1)

        # times when the Equation of time is maximum/minimum
        self.assertAlmostEqual(get_sun(Time("2020-02-11 12:14:15")).transform_to(I3Dir()).az.degree, 90, 1)
        self.assertAlmostEqual(get_sun(Time("2020-05-14 11:56:19")).transform_to(I3Dir()).az.degree, 90, 1)
        self.assertAlmostEqual(get_sun(Time("2020-07-26 12:06:36")).transform_to(I3Dir()).az.degree, 90, 1)
        self.assertAlmostEqual(get_sun(Time("2020-11-03 11:43:35")).transform_to(I3Dir()).az.degree, 90, 1)

        self.assertAlmostEqual(get_sun(Time("2020-03-20 03:50")).transform_to(I3Dir()).zen.degree, 90, 1)
        self.assertAlmostEqual(
            get_sun(Time("2020-06-20 21:43")).transform_to(I3Dir()).zen.degree, 113.44, 1
        )
        self.assertAlmostEqual(get_sun(Time("2020-09-22 13:31")).transform_to(I3Dir()).zen.degree, 90, 1)
        self.assertAlmostEqual(get_sun(Time("2020-12-21 10:03")).transform_to(I3Dir()).zen.degree, 66.56, 1)

    def test_sun_array(self):

        t = Time("2020-03-20 0:00")
        x1 = np.arange(366)
        t1 = t + x1 * day
        sun1 = get_sun(t1).transform_to(I3Dir())
        ref1 = I3Dir(zen=(90 + 23.44 * np.sin(x1 / len(x1) * 2 * np.pi)) * deg, az=270 * deg)
        assert_allclose(sun1.zen.degree, ref1.zen.degree, rtol=0.02)
        assert_allclose(sun1.az.degree, ref1.az.degree, rtol=0.02)
        assert_allclose(0, sun1.separation(ref1).degree, atol=4.2)

        x2 = np.linspace(0, 1, 1441)
        t2 = t + x2 * day
        sun2 = get_sun(t2).transform_to(I3Dir())
        ref2 = I3Dir(zen=90 * deg, az=(268.2 + x2 * 360) * deg)
        assert_allclose(sun2.zen.degree, ref2.zen.degree, rtol=0.004)
        assert_allclose((sun2.az - ref2.az).wrap_at(180 * deg).degree, 0, atol=0.1)
        assert_allclose(0, sun2.separation(ref2).degree, atol=0.4)

    @unittest.skipIf("astro" not in globals(), "Not in an icetray environment")
    def test_icetray(self):

        t0 = Time("2020-03-20 0:00")
        t = t0 + np.arange(366) * u.day

        np.random.seed(0)
        zen = 1 * np.pi * np.random.random(len(t))
        azi = 2 * np.pi * np.random.random(len(t))

        i3dir = I3Dir(zen=zen * u.radian, az=azi * u.radian, obstime=t)
        equa = i3dir.transform_to(ICRS())
        ra, dec = astro.dir_to_equa(zen, azi, t.mjd)

        assert_allclose(equa.ra.radian, ra, atol=1e-3)
        assert_allclose(equa.dec.radian, dec, atol=1e-5)

        crab = SkyCoord.from_name("Crab")
        i3crab = crab.transform_to(I3Dir(obstime=t))
        zenith, azimuth = astro.equa_to_dir(crab.ra.radian, crab.dec.radian, t.mjd)
        assert_allclose(zenith, i3crab.zen.radian, atol=1e-5)
        assert_allclose(azimuth, i3crab.az.radian, atol=2e-5)

        i3sun = get_sun(t).transform_to(I3Dir())
        sun_zen, sun_azi = astro.sun_dir(t.mjd)
        assert_allclose(sun_zen, i3sun.zen.radian, atol=2e-5)
        assert_allclose(sun_azi, i3sun.az.radian, atol=1e-4)

        i3moon = get_moon(t).transform_to(I3Dir())
        moon_zen, moon_azi = astro.moon_dir(t.mjd)
        assert_allclose(moon_zen, i3moon.zen.radian, atol=1e-4)
        assert_allclose(moon_azi, i3moon.az.radian, atol=1e-4)


if __name__ == "__main__":
    unittest.main()
