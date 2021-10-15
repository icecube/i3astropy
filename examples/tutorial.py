#!/usr/bin/env python

from astropy.coordinates import ICRS, get_moon
from astropy.time import Time
from astropy.units import radian

from i3astro import I3Dir

# Create an astropy time from a daq time
t = Time(2021, 155525688461058500, format="i3time")
print(t.iso)

# Convert IceCube coordinates to celestial
print(I3Dir(zen=0.73884 * radian, az=2.7348 * radian, obstime=t).transform_to(ICRS()))

# Convert Crab Neebula to IceCube coordinates
print(ICRS("05h34m31.94", "+22d00m52.2").transform_to(I3Dir(obstime=t)))

# Get the Moon location in IceCube Coordinates
print(get_moon(t).transform_to(I3Dir()))
