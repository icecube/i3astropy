#!/usr/bin/env python

# Copyright (c) 2022, The IceCube Collaboration, All Rights Reserved
# SPDX-License-Identifier:  GPL-3.0-or-later
# Author: Kevin Meagher

import warnings

import numpy as np
import pylab as plt
from astropy.coordinates import get_moon
from astropy.time import Time
from astropy.units import day
from matplotlib.dates import DateFormatter, YearLocator
from matplotlib.ticker import FormatStrFormatter

from i3astropy import I3Dir

warnings.simplefilter("ignore")

t0 = Time("2020-01-02 00:00")
T = 27.322
N_PERIODS = 252
toff = np.linspace(0, T, 50)

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=[12, 4.8])
min_t = []
min_zenith = []
min_dir = []

for i in range(N_PERIODS):
    t1 = t0 + i * T * day
    zenith = get_moon(t1 + toff * day).transform_to(I3Dir()).zen.degree
    # zenith = np.full(len(toff), 66)
    ax1.plot(toff, zenith)
    imin = np.argmin(zenith)
    min_t.append(t1 + toff[imin] * day)
    min_zenith.append(zenith[imin])
    print(f"{i:03}", min_t[-1], min_zenith[-1])

ax1.set_xlabel("Phase of Moon [Days]")
ax1.set_xlim(0, T)
ax1.set_ylim(120, 60)
ax1.set_ylabel("Moon Zenith")
ax1.yaxis.set_major_formatter(FormatStrFormatter("%d\u00b0"))
ax1.grid()

min_t = Time(min_t)
ax2.plot(min_t.plot_date, min_zenith)

ax2.xaxis.set_major_formatter(DateFormatter("%Y"))
ax2.xaxis.set_major_locator(YearLocator(2))
ax2.set_xlim(Time("2020-01-01 00:00").plot_date, min_t[-1].plot_date)
ax2.set_xlabel("Year")
ax2.set_ylim(74, 61)
ax2.set_ylabel("Moon Zenith at Minimum")
ax2.yaxis.set_major_formatter(FormatStrFormatter("%d\u00b0"))
ax2.grid()

fig1.tight_layout()
plt.show()
