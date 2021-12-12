[![Tests](https://github.com/icecube/i3astro/actions/workflows/unit_tests.yml/badge.svg)](https://github.com/icecube/i3astro/actions/workflows/unit_tests.yml)
[![codecov](https://codecov.io/gh/icecube/i3astro/branch/main/graph/badge.svg?token=VSU1VR44Y2)](https://codecov.io/gh/icecube/i3astro)

# i3astro

Astropy support for IceCube local coordinates.

Astropy is an extensive software collection which includes routines to convert between local and celestial coordinates.
Because of IceCube's unique location at the Geographic South Pole, it doesn't use the standard altitude/azimuth
coordinate system used by astropy and other astronomical software packages.
Instead IceCube uses a zenith/azimuth where the azimuth is defined with respect to polar grid directions rather than the local meridian.
IceCube's local coordinate system is described in depth [here](https://docs.icecube.aq/icetray/main/projects/dataclasses/coordinates.html)).
This python module provides a small extension to astropy with a class called `i3astro.I3Dir`, an instance of
`astropy.coordinates.baseframe.BaseCoordinateFrame`, which can be used for all astropy coordinate transforms.
In addition, a time format is provided to allow IceCube DAQ times to be entered directly.
For more information on astropy coordinate transform see the [astropy docs](https://docs.astropy.org/en/stable/coordinates/index.html).

## Installation

To install with pip:

    pip install [--user] git+https://github.com/icecube/i3astro.git

If you want to develop i3astro you can install directly with flit.
The ``-s`` option will symlink the module into site-packages rather than copying it,
so that you can test changes without reinstalling the module:

    git clone git@github.com:icecube/i3astro.git
    cd i3astro
    flit install [--user] -s


## Tutorial

This example takes an event which occurred at DAQ time (2021, 155525688461058500) and was reconstructed in the direction of zenith = 0.73884 rad, azimuth = 2.7348 rad and converts it to RA and Dec.
In addition, get the moon position in IceCube coordinates at the same time.

```python

>>> from astropy.coordinates import ICRS, get_sun
>>> from astropy.time import Time
>>> from astropy.units import radian
>>> from i3astro import I3Dir

>>> # Create an astropy time from a daq time
>>> t = Time(2021, 155525688461058500, format="i3time")
>>> t.iso
'2021-06-30 00:09:28.846'

>>> # Convert IceCube coordinates to celestial
>>> I3Dir(zen=0.73884 * radian, az=2.7348 * radian, obstime=t).transform_to(ICRS())
<ICRS Coordinate: (ra, dec) in deg
    (213.61898984, -47.55991453)>

>>> # Convert Crab Nebula to IceCube coordinates
>>> ICRS("05h34m31.94", "+22d00m52.2").transform_to(I3Dir(obstime=t))
<I3Dir Coordinate (obstime=(2021, 155525688461058496)): (az, zen, r) in (deg, deg, )
    (286.7095864, 112.03182914, 1.)>

>>> # Get the Moon location in IceCube Coordinates
>>> get_sun(t).transform_to(I3Dir())
<SkyCoord (I3Dir: obstime=(2021, 155525688461058496)): (az, zen, r) in (deg, deg, )
    (271.46065373, 113.17296322, 1.)>

```
