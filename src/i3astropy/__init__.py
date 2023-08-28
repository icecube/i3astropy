# SPDX-FileCopyrightText: © 2022 the i3astropy contributors
#
# SPDX-License-Identifier: BSD-2-Clause

"""Astropy support for the IceCube Coordinate System."""

__all__ = ["I3Time", "I3Dir", "I3DirToAltAz", "AltAzToI3Dir", "i3location"]
__version__ = "0.1"

import erfa
import numpy as np
from astropy.coordinates import AltAz, EarthLocation
from astropy.coordinates.attributes import TimeAttribute
from astropy.coordinates.baseframe import BaseCoordinateFrame, RepresentationMapping, frame_transform_graph
from astropy.coordinates.representation import PhysicsSphericalDifferential, PhysicsSphericalRepresentation
from astropy.coordinates.transformations import CoordinateTransform
from astropy.time import TimeFormat
from astropy.time.core import ScaleValueError
from astropy.units import deg, m

# This is the nominal location of the center of the IceCube detector
i3location = EarthLocation(lat=-89.9944 * deg, lon=-62.6081 * deg, height=883.9 * m)


class I3Time(TimeFormat):
    """Time format for IceCube's DAQ format.

    Times are expressed as two integers: The UTC year and the number of
    tenths of nanoseconds since the start of the UTC Year including leap seconds.

    Notes
    -----
        :py:class:`astropy.Time` internal storage format has higher precision than DAQ ticks.
        However, when initializing Time objects astropy converts all parameters to float64,
        which for values DAQ times close to the end of the year can result in a loss of
        precision of up to 64 DAQ ticks (6.4 nanoseconds).
    """

    name = "i3time"  # Unique format name
    _default_scale = "utc"
    _DAQ_TICKS_IN_DAY = int(864e12)

    def set_jds(self, val1, val2):
        """Set the internal jd1 and jd2 values from the input year and DAQ time.

        Year is in val1 and daq time is in val2.
        The input values are expected to conform to this format, as
        validated by self._check_val_type(val1, val2) during __init__.
        """
        self._scale = self._check_scale(self._scale)  # Validate scale.
        if self._scale != "utc":
            msg = f"Got scale '{self._scale}', The only allowed scale for I3Time is 'utc'"
            raise ScaleValueError(msg)

        year = val1.astype(int)
        if np.any(year != val1):
            msg = f"expected integer type for year, got {val1}"
            raise ValueError(msg)

        # divide DAQ time into days and remainder of days
        days, remainder = np.divmod(val2, self._DAQ_TICKS_IN_DAY)
        # convert start of year to tai
        tai1, tai2 = erfa.utctai(*erfa.dtf2d(b"UTC", year, 1, 1, 0, 0, 0))
        # add daq time to start of year in tai and then convert to utc
        self.jd1, self.jd2 = erfa.taiutc(tai1 + days, tai2 + remainder / self._DAQ_TICKS_IN_DAY)

    @property
    def value(self):
        """Return year and daq time from internal jd1, jd2."""
        assert self.scale == "utc"
        # create empty output array with correct type
        out = np.empty(self.jd1.shape, dtype=[("year", "i4"), ("daq_time", "i8")])
        # get the current year
        out["year"], _, _, _ = erfa.d2dtf(b"UTC", 10, self.jd1, self.jd2_filled)
        # get the start of the year in tai jd
        y1, y2 = erfa.utctai(*erfa.dtf2d(b"UTC", out["year"], 1, 1, 0, 0, 0))
        # get the current time in tai
        t1, t2 = erfa.utctai(self.jd1, self.jd2)
        # calculate time deltas and convert to integer
        d1 = np.round((t1 - y1) * self._DAQ_TICKS_IN_DAY).astype(np.int64)
        d2 = np.round((t2 - y2) * self._DAQ_TICKS_IN_DAY).astype(np.int64)
        # add time deltas
        out["daq_time"] = d1 + d2
        # return masked array
        return self.mask_if_needed(out.view(np.recarray))


class I3Dir(BaseCoordinateFrame):
    """Representation the standard IceCube direction coordinates: zenith and azimuth angles.

    Zenith angle is the angle between the direction of the particle and the vertical (+z) direction.
    Azimuth angle is the projection in the xy-plane measured counterclockwise from the +x direction
    such that the +y direction is at azimuth=90°.
    """

    frame_specific_representation_info = {
        PhysicsSphericalRepresentation: [
            RepresentationMapping("theta", "zen"),
            RepresentationMapping("phi", "az"),
        ],
    }
    default_representation = PhysicsSphericalRepresentation
    default_differential = PhysicsSphericalDifferential

    obstime = TimeAttribute(default=None)
    location = i3location

    def __init__(self, *args, **kwargs) -> None:
        """Initialize I3Dir."""
        if "zen" in kwargs and "az" in kwargs and "r" not in kwargs:
            kwargs["r"] = 1
        super().__init__(*args, **kwargs)


class I3DirToAltAz(CoordinateTransform):
    """Convert from icecube coordinates to standard astronomy coordinates.

    For local directions.
    """

    def __init__(self) -> None:
        """Initialize I3DirToAltAz."""
        super().__init__(I3Dir, AltAz)

    def __call__(self, i3dir, *args, **kwargs):  # noqa: ARG002
        """Call I3DirToAltAz."""
        alt = 90 * deg - i3dir.zen
        azi = 90 * deg - i3dir.az - i3location.lon
        return AltAz(alt=alt, az=azi, location=i3dir.location, obstime=i3dir.obstime)


class AltAzToI3Dir(CoordinateTransform):
    """Convert from local astronomy coordinates to IceCube coordinates."""

    def __init__(self) -> None:
        """Initialize AltAzToI3Dir."""
        super().__init__(AltAz, I3Dir)

    def __call__(self, altaz, *args, **kwargs):  # noqa: ARG002
        """Call AltAzToI3Dir."""
        zen = 90 * deg - altaz.alt
        azi = 90 * deg - altaz.az - i3location.lon
        return I3Dir(az=azi, zen=zen, obstime=altaz.obstime)


frame_transform_graph.add_transform(I3Dir, AltAz, I3DirToAltAz())
frame_transform_graph.add_transform(AltAz, I3Dir, AltAzToI3Dir())
