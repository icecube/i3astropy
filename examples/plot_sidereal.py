#!/usr/bin/env python

# SPDX-FileCopyrightText: © 2026 the i3astropy contributors
#
# SPDX-License-Identifier: BSD-2-Clause

"""Plots example usage of sidereal functions for three months and the Fourier Transform."""

import numpy as np
import pylab as plt
from astropy.time import Time, TimeDelta
from matplotlib import dates as mdates
from matplotlib import ticker

from i3astropy import gmast, gmest, gmst, gmt

_DAYS_IN_YEAR = 365.2425
funcs = [
    (gmast, "Anti Sidereal"),
    (gmt, "Solar"),
    (gmst, "Sidereal"),
    (gmest, "Extended Sidereal"),
]

fig, axis = plt.subplots(2, 2, figsize=(9.6 * 2, 4.8 * 2))

for i, month in enumerate((8, 9, 10)):
    reftime = Time(f"2022-{month}-21 00:00:00", scale="utc")
    times = reftime + TimeDelta(np.linspace(-1.98, 2.98, 1000), format="jd")
    for func, name in funcs:
        y = func(times)
        print(type(y), y.unit)
        axis[i // 2, i % 2].plot(times.plot_date, y.hourangle, label=name)

    locator = mdates.DayLocator()
    axis[i // 2, i % 2].xaxis.set_major_locator(locator)
    axis[i // 2, i % 2].xaxis.set_minor_locator(mdates.HourLocator(range(0, 23, 6)))
    axis[i // 2, i % 2].xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
    axis[i // 2, i % 2].set_ylim(0, 24)
    axis[i // 2, i % 2].set_title(reftime.strftime("%B"))

DURATION = 36520 // 2
SAMPLES_PER_DAY = 100
SAMPLES = int(DURATION * SAMPLES_PER_DAY)

reftime = Time("1970-1-1 00:00:00", scale="utc")
times = reftime + TimeDelta(np.linspace(0, DURATION, SAMPLES), format="jd")
w = np.linspace(0, SAMPLES_PER_DAY // 2, SAMPLES // 2) * _DAYS_IN_YEAR

for func, name in funcs:
    y = func(times).hour
    f = np.abs(np.fft.rfft(y, norm="ortho"))[:-1]
    f[0] = 0
    axis[1, 1].plot(w, f, label=name)
    print(f"{name:>17} {w[f.argmax()]:10.6f} {w[1] - w[0]:8.6f}")
axis[1, 1].xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
axis[1, 1].legend(loc=1)
axis[1, 1].set_xlim(363.0, 368.5)
axis[1, 1].set_title("FFT")

fig.tight_layout()
fig.savefig("plot_sidereal.svg")
plt.show()
