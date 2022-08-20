#!/usr/bin/env python

# Copyright (c) 2022, The IceCube Collaboration, All Rights Reserved
# SPDX-License-Identifier:  GPL-3.0-or-later
# Author: Kevin Meagher

"""
takes a file from the Neutrino Sources Working Group's datasets and calculate
the Right Ascension and Declination of the events from the Zenith and Azimuth
and compare to the coordinates in the data file. It can be seen that i3astropy
agrees with previous methods to within 0.0002°, the high values for RA and
azimuth are caused by events close to the pole.
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import ICRS
from astropy.time import Time
from numpy import allclose

from i3astropy import I3Dir

# need file from /data/ana/analyses/ps_tracks/version-004-p00/IC86_2019_exp.npy
f = np.load("IC86_2019_exp.npy")
N = 10000

time = Time(f["time"][:N], scale="ut1", format="mjd")
i3dir = I3Dir(zen=f["zen"][:N] * u.rad, az=f["azi"][:N] * u.rad, obstime=time)
icrs = ICRS(ra=f["ra"][:N] * u.rad, dec=f["dec"][:N] * u.rad)
icrs2 = i3dir.transform_to(ICRS())
i3dir2 = icrs.transform_to(I3Dir(obstime=time))

assert allclose(f["ra"][:N], icrs2.ra.radian, atol=1e-3)
assert allclose(f["dec"][:N], icrs2.dec.radian, atol=1e-5)
assert allclose(icrs2.separation(icrs).radian, 0, atol=1e-5)
assert allclose(f["zen"][:N], i3dir2.zen.radian, atol=1e-5)
assert allclose(f["azi"][:N], i3dir2.az.radian, atol=2e-4)
assert allclose(i3dir.separation(i3dir2).radian, 0, atol=1e-5)

print(f"Checked {N} events")
print(f"Max diff RA            : {np.rad2deg(max(abs(f['ra'][:N] - icrs2.ra.radian))):4.6f}°")
print(f"Max diff dec           : {np.rad2deg(max(abs(f['dec'][:N] - icrs2.dec.radian))):4.6f}°")
print(f"Max diff sky separation: {max(icrs2.separation(icrs).degree):4.6f}°")
print(f"Max diff zenith        : {np.rad2deg(max(abs(f['zen'][:N] - i3dir2.zen.radian))):4.6f}°")
print(f"Max diff azimuth       : {np.rad2deg(max(abs(f['azi'][:N] - i3dir2.az.radian))):4.6f}°")
print(f"Max diff I3 separation : {max(i3dir.separation(i3dir2).degree):4.6f}°")
