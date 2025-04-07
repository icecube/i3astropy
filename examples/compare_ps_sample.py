#!/usr/bin/env python

# SPDX-FileCopyrightText: © 2022 the i3astropy contributors
#
# SPDX-License-Identifier: BSD-2-Clause

"""Compares the RA and Dec from a Nu Sources Dataset with the output if i3astropy.

Takes a file from the Neutrino Sources Working Group's datasets and calculate
the Right Ascension and Declination of the events from the Zenith and Azimuth
in that file and compares it to the RA and Dec in the data file.
It can be seen that i3astropy agrees with previous methods to within 0.0002°,
the high values for RA and azimuth are caused by events close to the pole.
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

ic_ra = np.cos(f["dec"][:N]) * f["ra"][:N]
ic_az = np.sin(f["zen"][:N]) * f["azi"][:N]
ap_ra = np.cos(icrs2.dec.radian) * icrs2.ra.radian
ap_az = np.sin(i3dir2.zen.radian) * i3dir2.az.radian

assert allclose(ic_ra, ap_ra, atol=1e-5)
assert allclose(f["dec"][:N], icrs2.dec.radian, atol=1e-5)
assert allclose(icrs2.separation(icrs).radian, 0, atol=1e-5)
assert allclose(f["zen"][:N], i3dir2.zen.radian, atol=1e-5)
assert allclose(ic_az, ap_az, atol=1e-5)
assert allclose(i3dir.separation(i3dir2).radian, 0, atol=1e-5)

print(f"Checked {N} events")
print(f"Max diff RA            : {np.rad2deg(max(abs(ic_ra - ap_ra))):4.6f}°")
print(f"Max diff dec           : {np.rad2deg(max(abs(f['dec'][:N] - icrs2.dec.radian))):4.6f}°")
print(f"Max diff sky separation: {max(icrs2.separation(icrs).degree):4.6f}°")
print(f"Max diff zenith        : {np.rad2deg(max(abs(f['zen'][:N] - i3dir2.zen.radian))):4.6f}°")
print(f"Max diff azimuth       : {np.rad2deg(max(abs(ic_az - ap_az))):4.6f}°")
print(f"Max diff I3 separation : {max(i3dir.separation(i3dir2).degree):4.6f}°")
