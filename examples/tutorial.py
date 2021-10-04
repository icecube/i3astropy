#!/usr/bin/env python

from astropy.units import radian
from astropy.coordinates import ICRS, get_moon
from astropy.time import Time
from i3astro import I3Dir

t= Time(2021,155525688461058500,format='i3time')
print(t.iso)
print(I3Dir(zen=0.73884*radian, az=2.7348*radian,obstime=t).transform_to(ICRS()))
print (get_moon(t).transform_to(I3Dir()))
