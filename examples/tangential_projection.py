# SPDX-FileCopyrightText: 2024 the i3astropy contributors
#
# SPDX-License-Identifier: BSD-2-Clause

"""example to show how to make plots in a tangential projection using astropy's wcs."""

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.wcs import WCS
from matplotlib.colors import LogNorm

from i3astropy import I3Dir

# this is were we will center our plot
center = SkyCoord(90 * u.deg, 80 * u.deg)
# the is the size across the plot in degree
SIZE = 10

wcs_helix_list = WCS(naxis=2)
# for some reason the code is 1 = bottom left, 2= bottom right 1.5 = center
wcs_helix_list.wcs.crpix = [1.5, 1.5]
# set the center of the plot
wcs_helix_list.wcs.crval = [center.ra.deg, center.dec.deg]
# use degrees
wcs_helix_list.wcs.cunit = ["deg", "deg"]
wcs_helix_list.wcs.ctype = ["RA---TAN", "DEC--TAN"]
# make the plot approx 10 degree across
wcs_helix_list.wcs.cdelt = [SIZE, SIZE]

# load the datafile but cut it to 30000 events so we can see what is happening
f = np.load("IC86_2019_exp.npy")[:30000]

# create a figure
fig = plt.figure(figsize=(12, 10))
ax = plt.subplot(projection=wcs_helix_list)

# save the limits of the plot for later
xlim = ax.get_xlim()
ylim = ax.get_ylim()

# set the labels and make a grid
ax.set_xlabel("RA")
ax.set_ylabel("Dec")
ax.grid(color="grey", ls="dotted")


time = Time(f["time"], scale="ut1", format="mjd")
i3dir = I3Dir(zen=f["zen"] * u.rad, az=f["azi"] * u.rad, obstime=time)
icrs = SkyCoord(ra=f["ra"] * u.rad, dec=f["dec"] * u.rad)

MAX_ANG_ERR = 2
# Cut out events which are too far from the center of our plot or have big errors
cut = (icrs.separation(center) < SIZE / 1.5 * u.deg) * (np.rad2deg(f["angErr"]) < MAX_ANG_ERR)

# make a color scale based on energy
cmap = plt.cm.rainbow
norm = LogNorm(vmin=1e2, vmax=1e4)
colors = cmap(norm(10 ** (f["logE"][cut])))

# plot our events based using coordinate transform from zenith-azimuth with a plut sign
icrs2 = SkyCoord(i3dir[cut], frame="i3dir").transform_to("icrs")
ax.scatter(
    icrs2.ra.deg, icrs2.dec.deg, transform=ax.get_transform("icrs"), s=100, facecolor=colors, marker="+"
)

# plot our events based on the RA-Dec in the file with a cross overlaying the plus signs we just made
icrs = icrs[cut]
ax.scatter(
    icrs.ra.deg, icrs.dec.deg, transform=ax.get_transform("icrs"), s=100, facecolor=colors, marker="x"
)

# Plot a Circle around each event with radius of its angErr
for i, coord in enumerate(icrs):
    ax.add_patch(
        SphericalCircle(
            SkyCoord(coord),
            f["angErr"][cut][i] * u.rad,
            edgecolor=colors[i],
            facecolor="none",
            transform=ax.get_transform("icrs"),
        )
    )

# set the limits back to what they were before we plotted stuff off screen
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)

# add color bar to the side of the plot
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cbar = fig.colorbar(sm, ax=ax)
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel("Energy [GeV]", rotation=270)

fig.tight_layout()
plt.show()
