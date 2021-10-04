#!/usr/bin/env python
# This example can verify that the azimuth is defined correctly with reference to the sun
# It is obvious that the sun should be at Grid South at midnight,
# Grid East at 6:00, Grid North at 12:00, and Grid West at 18:00
# the example makes a plot to verify that coordinates are defined correctly

import matplotlib.dates as mdates
import numpy as np
import pylab as plt
from astropy.coordinates import ICRS, get_sun
from astropy.time import Time
from astropy.units import day, deg
from matplotlib import patches, ticker

from i3astro import I3Dir


def draw_circ(ax, radius, centX, centY, angle_, theta2_, color_="black"):
    arc = patches.Arc(
        [centX, centY],
        radius,
        radius,
        angle=angle_,
        theta1=0,
        theta2=theta2_,
        capstyle="round",
        linestyle="-",
        lw=1.5,
        color=color_,
    )
    ax.add_patch(arc)

    endX = centX + (radius / 2) * np.cos(np.deg2rad(theta2_ + angle_))  # Do trig to determine end position
    endY = centY + (radius / 2) * np.sin(np.deg2rad(theta2_ + angle_))

    ax.add_patch(  # Create triangle as arrow head
        patches.RegularPolygon(
            (endX, endY),  # (x,y)
            3,  # number of vertices
            radius / 7,  # radius
            np.deg2rad(angle_ + theta2_),  # orientation
            color=color_,
        )
    )


fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=[12, 4.8], gridspec_kw=dict(width_ratios=(1, 1.1), wspace=0.35, left=0.06, right=0.97)
)

acolor = "tab:blue"

draw_circ(ax1, 200, 0, 0, 0, 315, acolor)
ax1.text(
    100 / 2 ** 0.5,
    -100 / 2 ** 0.5,
    "Azimuth",
    horizontalalignment="left",
    verticalalignment="top",
    color=acolor,
    size=15,
)
ax1.plot([-200, 200], [0, 0], color=acolor)
ax1.plot([0, 0], [-200, 200], color=acolor)

xtoe = lambda x: x / 0.3048 + 46500
etox = lambda E: 0.3048 * (E - 46500)
yton = lambda y: y / 0.3048 + 52200
ntoy = lambda N: 0.3048 * (N - 52200)


S, E, N, _ = np.loadtxt("IceCubeAsBuiltHoleCoordinates.txt").T

X = []
Y = []
factor = 1.14
for string in [1, 6, 50, 74, 72, 78, 75, 31, 1]:
    i = np.argwhere(S == string)[0][0]
    X.append(factor * etox(E[i]))
    Y.append(factor * ntoy(N[i]))

ax1.plot(X, Y, color="k")
ax1.scatter(etox(E), ntoy(N), color="k", marker=".", s=1.5)

secxax = ax1.secondary_xaxis("top", functions=(xtoe, etox))
secxax.xaxis.set_major_locator(ticker.MultipleLocator(1000))
secxax.xaxis.set_major_formatter("E {x:0,.0f}′")
secyax = ax1.secondary_yaxis("right", functions=(yton, ntoy))
secyax.yaxis.set_major_locator(ticker.MultipleLocator(1000))
secyax.yaxis.set_major_formatter("N {x:0,.0f}′")
ax1.set_aspect("equal")

ax1.set_xlim(-900, 900)
ax1.set_ylim(-700, 700)
ax1.set_xlabel("x [m]")
ax1.set_ylabel("y [m]")

for x, y, t in [(0, -630, 0), (750, 0, 6), (0, 630, 12), (-750, 0, 18)]:
    ax1.text(
        x,
        y,
        f"\u2600{t:02d}:00  ",
        c="tab:orange",
        fontsize=12,
        horizontalalignment="center",
        verticalalignment="center",
    )

t = Time("2021-03-20") + np.linspace(0, 1, 1001) * day
fpa_azimuth = ICRS(ra=0 * deg, dec=0 * deg).transform_to(I3Dir(obstime=t)).az
sun_azimuth = get_sun(t).transform_to(I3Dir()).az

plt.plot(t.plot_date, np.rad2deg(fpa_azimuth), ls="-", label="First Point of Ares")
plt.plot(t.plot_date, np.rad2deg(sun_azimuth), ls="--", label="Sun")

ax2.xaxis.set_major_locator(mdates.AutoDateLocator())
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax2.xaxis.set_minor_locator(mdates.HourLocator())

ax2.yaxis.set_major_locator(ticker.MultipleLocator(45))
ax2.yaxis.set_major_formatter("{x:0.0f}°")
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(15))

secax = ax2.secondary_yaxis("right")
secax.yaxis.set_major_locator(ticker.MultipleLocator(90))
secax.yaxis.set_major_formatter(
    ticker.FuncFormatter(lambda x, y: ["E", "N", "W", "S"][(int(x) % 360) // 90])
)
ax2.set_ylabel("Azimuth")
ax2.set_xlabel("UTC Time")
ax2.set_xlim(t[0].plot_date, t[-1].plot_date)
ax2.set_ylim(0, 360)

ax2.set_title("Vernal Equinox " + t[0].to_value("iso", subfmt="date"))
ax2.legend()
ax2.grid()
plt.show()
