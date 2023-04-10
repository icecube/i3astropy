#!/usr/bin/env python

# SPDX-FileCopyrightText: © 2022 the i3astropy contributors
#
# SPDX-License-Identifier: BSD-2-Clause

"""plots the IceCube Azimuth of the Sun and First Point of Ares on the Vernal
Equinox next to a diagram of the IceCube detector with IceCube Azimuth labeled
and the positions of the sun at certain times. This Demonstrates that
i3astropy is correctly orienting the IceCube coordinate system.
"""
import matplotlib.dates as mdates
import numpy as np
import pylab as plt
from astropy.coordinates import get_sun
from astropy.time import Time
from astropy.units import day
from i3astropy import I3Dir
from matplotlib import ticker

LEFT_MARGIN, RIGHT_MARGIN = 0.08, 0.02
ANALAMA_WIDTH = 0.2
TOP_MARGIN, BOTTOM_MARGIN = 0.01, 0.08
HORIZOTAL_GAP, VERTICAL_GAP = 0.1, 0

rect_analemma = [LEFT_MARGIN, BOTTOM_MARGIN, ANALAMA_WIDTH, 1 - BOTTOM_MARGIN - TOP_MARGIN]
rect_elevation = [
    LEFT_MARGIN + ANALAMA_WIDTH + HORIZOTAL_GAP,
    (1 + BOTTOM_MARGIN - TOP_MARGIN + HORIZOTAL_GAP) / 2,
    1 - LEFT_MARGIN - ANALAMA_WIDTH - HORIZOTAL_GAP - RIGHT_MARGIN,
    (1 - BOTTOM_MARGIN - TOP_MARGIN - HORIZOTAL_GAP) / 2,
]
rect_time = [
    LEFT_MARGIN + ANALAMA_WIDTH + HORIZOTAL_GAP,
    BOTTOM_MARGIN,
    1 - LEFT_MARGIN - ANALAMA_WIDTH - HORIZOTAL_GAP - RIGHT_MARGIN,
    (1 - BOTTOM_MARGIN - TOP_MARGIN - HORIZOTAL_GAP) / 2,
]

fig = plt.figure(figsize=(10, 6))
ax1 = fig.add_axes(rect_analemma)
ax2 = fig.add_axes(rect_elevation)
ax3 = fig.add_axes(rect_time, sharex=ax2)

T0 = Time("2021-01-01 12:00:00", format="iso", scale="utc")
T = T0 + np.arange(0, 366) * day
sun = get_sun(T).transform_to(I3Dir())

ax2.plot(T.plot_date, sun.zen.degree, ls="-", marker=None)
ax2.set_ylabel("Zenith")
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d\u00b0"))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(2))
ax2.set_ylim(116, 64)
ax2.grid()

ax3.plot(T.plot_date, sun.az.degree, ls="-", marker=None)
ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d\u00b0"))
ax3.yaxis.set_minor_locator(ticker.MultipleLocator(1))
ax3.set_ylabel("Azimuth")
ax3.set_ylim(86, 95)
ax3.grid()

ax3.xaxis.set_major_locator(mdates.MonthLocator([1, 3, 5, 7, 9, 11]))
ax3.xaxis.set_major_formatter(mdates.DateFormatter("1 %b"))
ax3.xaxis.set_minor_locator(mdates.MonthLocator())
ax3.set_xlim((T.plot_date[0] - 1, T.plot_date[-1]))

plt_info = ax1.plot(sun.az, sun.zen)
for n, d in enumerate(T[:-1]):
    if d.datetime.day == 1:
        slope = -(sun[n + 1].az.degree - sun[n - 1].az.degree) / (
            sun[n + 1].zen.degree - sun[n - 1].zen.degree
        )
        SIZE = 1
        x1 = -SIZE / (1 + slope**2) ** 0.5
        x2 = +SIZE / (1 + slope**2) ** 0.5

        ax1.plot(
            [sun[n].az.degree + x1, sun[n].az.degree + x2],
            [
                sun[n].zen.degree + slope * x1,
                sun[n].zen.degree + slope * x2,
            ],
            c=plt_info[0].get_c(),
        )

        if sun[n].az.degree < 90:
            HALIGN = "left"
            text_az = sun[n].az.degree - 1
            x = -1.5 * SIZE / (1 + slope**2) ** 0.5
            txt_az = sun[n].az.degree + x
            txt_zen = sun[n].zen.degree + slope * x
        else:
            HALIGN = "right"
            text_az = sun[n].az.degree + 1
            x = 1.5 * SIZE / (1 + slope**2) ** 0.5
            txt_az = sun[n].az.degree + x
            txt_zen = sun[n].zen.degree + slope * x
        ax1.text(
            txt_az,
            txt_zen,
            "1 " + d.strftime("%b"),
            verticalalignment="center",
            horizontalalignment=HALIGN,
        )

ax1.set_xlabel("Azimuth")
ax1.set_xlim(100, 80)
ax1.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d\u00b0"))
ax1.set_ylabel("Zenith")
ax1.set_ylim(116, 64)
ax1.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(1))
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d\u00b0"))

plt.show()
