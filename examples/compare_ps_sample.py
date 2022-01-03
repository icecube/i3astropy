#!/usr/bin/env python
import numpy as np
from astropy import units as u
from astropy.coordinates import ICRS
from astropy.time import Time
from numpy import allclose

from i3astropy import I3Dir

#need file from /data/ana/analyses/ps_tracks/version-004-p00/IC86_2019_exp.npy
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
print("Max diff RA            : {:4.6f}°".format(np.rad2deg(max(abs(f["ra"][:N] - icrs2.ra.radian)))))
print("Max diff dec           : {:4.6f}°".format(np.rad2deg(max(abs(f["dec"][:N] - icrs2.dec.radian)))))
print("Max diff sky separation: {:4.6f}°".format(max(icrs2.separation(icrs).degree)))
print("Max diff zenith        : {:4.6f}°".format(np.rad2deg(max(abs(f["zen"][:N] - i3dir2.zen.radian)))))
print("Max diff azimuth       : {:4.6f}°".format(np.rad2deg(max(abs(f["azi"][:N] - i3dir2.az.radian)))))
print("Max diff I3 separation : {:4.6f}°".format(max(i3dir.separation(i3dir2).degree)))
