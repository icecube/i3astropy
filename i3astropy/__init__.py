"""
Astropy support for IceCube local coordinates.
"""

__all__ = ["I3Time", "I3Dir", "i3location"]
__version__ = "0.1"

import erfa
import numpy as np
from astropy.coordinates import AltAz, EarthLocation
from astropy.coordinates.attributes import TimeAttribute
from astropy.coordinates.baseframe import BaseCoordinateFrame, RepresentationMapping, frame_transform_graph
from astropy.coordinates.representation import PhysicsSphericalDifferential, PhysicsSphericalRepresentation
from astropy.coordinates.transformations import CoordinateTransform
from astropy.time import Time, TimeFormat
from astropy.time.core import ScaleValueError
from astropy.units import deg, m, picosecond, second

# This is the nominal location of the center of the IceCube detector
i3location = EarthLocation(lat=-89.9944 * deg, lon=-62.6081 * deg, height=883.9 * m)


class I3Time(TimeFormat):
    """
    Time format for IceCube's DAQ format, times are expressed as two integers:
    The UTC year and the number of tenths of nanoseconds since the start of the UTC Year
    including leapseconds
    """

    name = "i3time"  # Unique format name
    _default_scale = "utc"
    _DAQ_SEC = int(1e10)

    def set_jds(self, val1, val2):
        """
        Set the internal jd1 and jd2 values from the input val1, val2.
        The input values are expected to conform to this format, as
        validated by self._check_val_type(val1, val2) during __init__.
        """

        self._scale = self._check_scale(self._scale)  # Validate scale.
        if self._scale != "utc":
            raise ScaleValueError(f"Got scale '{self._scale}', The only allowed scale for I3Time is 'utc'")

        year = val1.astype(int)
        if np.any(year != val1):
            raise ValueError(f"expected integer type for year, got {val1}")
        sec, daq = np.divmod(val2, self._DAQ_SEC)
        time = Time({"year": year}, format="ymdhms", scale="utc")
        time += sec * second + daq * 100 * picosecond
        self.jd1, self.jd2 = (time.jd1, time.jd2)

    @property
    def value(self):
        """
        Return format ``value`` property from internal jd1, jd2
        """

        out = np.empty(self.jd1.shape, dtype=[("year", "i4"), ("daq_time", "i8")])
        scale = self.scale.upper().encode("ascii")
        out["year"], _, _, _ = erfa.d2dtf(scale, 10, self.jd1, self.jd2_filled)
        time = Time(self.jd1, self.jd2, scale="utc", format="jd")
        year_start = Time({"year": out["year"]}, format="ymdhms", scale="utc")
        out["daq_time"] = (time - year_start).to_value(100 * picosecond).round().astype(int)
        out = out.view(np.recarray)
        return self.mask_if_needed(out)


class I3Dir(BaseCoordinateFrame):
    """
    Representation the standard IceCube direction coordinates: zenith and azimuth angles

    Zenith angle is the angle between the direction of the particle and the vertical (+z) direction.
    Azimuth angle is the projection in the xy-plane measured counterclockwise from the +x direction
    such that the +y direction is at azimuth=90°.
    """

    frame_specific_representation_info = {
        PhysicsSphericalRepresentation: [
            RepresentationMapping("theta", "zen"),
            RepresentationMapping("phi", "az"),
        ]
    }
    default_representation = PhysicsSphericalRepresentation
    default_differential = PhysicsSphericalDifferential

    obstime = TimeAttribute(default=None)
    location = i3location

    def __init__(self, *args, **kwargs):
        if "zen" in kwargs and "az" in kwargs and "r" not in kwargs:
            kwargs["r"] = 1
        super().__init__(*args, **kwargs)


class I3DirToAltAz(CoordinateTransform):
    """
    Convert from icecube coordinates to standard astronomy coordinates for
    local directions
    """

    def __init__(self):
        super().__init__(I3Dir, AltAz)

    def __call__(self, i3dir, *args, **kwargs):

        alt = 90 * deg - i3dir.zen
        azi = 90 * deg - i3dir.az - i3location.lon
        return AltAz(alt=alt, az=azi, location=i3dir.location, obstime=i3dir.obstime)


class AltAzToI3Dir(CoordinateTransform):
    """
    Convert from local astronomy coordinates to IceCube coordinates
    """

    def __init__(self):
        super().__init__(AltAz, I3Dir)

    def __call__(self, altaz, *args, **kwargs):
        zen = 90 * deg - altaz.alt
        azi = 90 * deg - altaz.az - i3location.lon
        return I3Dir(az=azi, zen=zen, obstime=altaz.obstime)


frame_transform_graph.add_transform(I3Dir, AltAz, I3DirToAltAz())
frame_transform_graph.add_transform(AltAz, I3Dir, AltAzToI3Dir())
