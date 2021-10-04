#!/usr/bin/env python

import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import numpy as np
import pylab as plt
from astropy.coordinates import get_sun
from astropy.time import Time
from astropy.units import day

from i3astro import I3Dir

lmargin, rmargin = 0.08, 0.02
awidth = 0.2
tmargin, bmargin = 0.01, 0.08
hgap, vgap = 0.1, 0

rect_analemma = [lmargin, bmargin, awidth, 1 - bmargin - tmargin]
rect_elevation = [
    lmargin + awidth + hgap,
    (1 + bmargin - tmargin + hgap) / 2,
    1 - lmargin - awidth - hgap - rmargin,
    (1 - bmargin - tmargin - hgap) / 2,
]
rect_time = [
    lmargin + awidth + hgap,
    bmargin,
    1 - lmargin - awidth - hgap - rmargin,
    (1 - bmargin - tmargin - hgap) / 2,
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
        size = 1
        x1 = -size / (1 + slope ** 2) ** 0.5
        x2 = +size / (1 + slope ** 2) ** 0.5

        ax1.plot(
            [sun[n].az.degree + x1, sun[n].az.degree + x2],
            [
                sun[n].zen.degree + slope * x1,
                sun[n].zen.degree + slope * x2,
            ],
            c=plt_info[0].get_c(),
        )

        if sun[n].az.degree < 90:
            halign = "left"
            text_az = sun[n].az.degree - 1
            x = -1.5 * size / (1 + slope ** 2) ** 0.5
            txt_az = sun[n].az.degree + x
            txt_zen = sun[n].zen.degree + slope * x
        else:
            halign = "right"
            text_az = sun[n].az.degree + 1
            x = 1.5 * size / (1 + slope ** 2) ** 0.5
            txt_az = sun[n].az.degree + x
            txt_zen = sun[n].zen.degree + slope * x
        ax1.text(
            txt_az, txt_zen, "1 " + d.strftime("%b"), verticalalignment="center", horizontalalignment=halign
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
