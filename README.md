# i3astro

An interface between astropy and IceCube's local coordinate system.

## Tutorial

This for an event which occured at daq time (2021, 155525688461058500) and was reconstructed in the direction of (0.73884, 2.7348) radians convert to RA and Dec.
In addition, get the moon position in IceCube coordinates at the same time.

    from astropy.units import radian
    from astropy.coordinates import ICRS, get_moon
    from astropy.time import Time
    from i3astro import I3Dir

    t= Time(2021,155525688461058500,format='i3time')
    print(t.iso)
    print(I3Dir(zen=0.73884*radian, az=2.7348*radian,obstime=t).transform_to(ICRS()))
    print (get_moon(t).transform_to(I3Dir()))
