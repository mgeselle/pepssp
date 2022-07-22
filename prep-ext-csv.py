#!/usr/bin/python3

from astropy.coordinates import SkyCoord, Longitude, Latitude, EarthLocation, AltAz
import astropy.time as atime
from csv import DictReader
from math import sqrt, log10, log, fabs
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
from statistics import mean, stdev
import sys


class Measurement:
    def __init__(self, star, time_start):
        self.star = star
        self.time = atime.Time(time_start, scale='utc')
        self.cnt = dict()
        self.std_cnt = dict()
        self.end_time = None
        self.airmass = None

    def add_counts(self, filte, integ_time, c1, c2, c3):
        divisor = integ_time / 1000
        net_cnts = [c1 / divisor, c2 / divisor, c3 / divisor]
        self.cnt[filt] = mean(net_cnts)
        self.std_cnt[filt] = stdev(net_cnts)

    def add_sky_counts(self, sky_time, filt, integ_time, c1, c2, c3):
        divisor = integ_time / 1000
        net_cnts = [c1 / divisor, c2 / divisor, c3 / divisor]
        sky_mean = mean(net_cnts)
        sky_stdev = stdev(net_cnts)
        self.cnt[filt] = self.cnt[filt] - sky_mean
        self.std_cnt[filt] = sqrt((self.std_cnt[filt]**2 + sky_stdev**2) / 2)
        self.end_time = atime.Time(sky_time, scale='utc')

    def calc_time(self):
        self.time = self.time + ((self.end_time - self.time) / 2)


coords = {
    'SAO 65890': SkyCoord('17h15m02.8343634s', '+36d48m32.98434s', frame='icrs'),
    'SAO 65921': SkyCoord('17h17m40.253198052s', '+37d17m29.422635288s', frame='icrs')
}


location = EarthLocation.from_geodetic(Longitude('8d32m12s'), Latitude('50d6m17s'), height=110.0)


measurements = []
cmp_star = None
pgm_star = None
num_pgm = 0
with open('/home/mgeselle/astro/pep-runs/20220718/second-order-ext.csv') as input_file:
    reader = DictReader(input_file)
    last_star = None
    current = None
    for row in reader:
        filt = row['Filter']
        star = row['StarId']
        if star != last_star:
            if current is not None:
                current.calc_time()
                measurements.append(current)
                if pgm_star is None:
                    pgm_star = star
            elif cmp_star is None:
                cmp_star = star
            current = Measurement(star, row['Timestamp'])
            last_star = star
            if star == pgm_star:
                num_pgm = num_pgm + 1
        if row['IsStar'] == 'True':
            current.add_counts(filt, int(row['IntegrationTime']),
                               int(row['Count1']), int(row['Count2']), int(row['Count3']))
        else:
            current.add_sky_counts(row['Timestamp'], filt, int(row['IntegrationTime']),
                                   int(row['Count1']), int(row['Count2']), int(row['Count3']))
    current.calc_time()
    measurements.append(current)

pre = None
mid = None
post = None

mag_cmp = dict()
mag_pgm = dict()
stdev_mag_cmp = dict()
stdev_mag_pgm = dict()
airmass = np.empty(num_pgm)
idx = 0
for measurement in measurements:
    if pre is not None and mid is None and measurement.star == cmp_star:
        pre = measurement
        continue
    if pre is None:
        pre = measurement
        continue
    elif mid is None:
        mid = measurement
        continue
    else:
        post = measurement
    mid_time = pre.time + (post.time - pre.time) / 2
    pos_cmp = coords[pre.star]
    altaz_cmp = pos_cmp.transform_to(AltAz(obstime=mid_time, location=location))
    pos_pgm = coords[mid.star]
    altaz_pgm = pos_pgm.transform_to(AltAz(obstime=mid_time, location=location))
    airmass[idx] = (altaz_cmp.secz + altaz_pgm.secz) / 2
    cnt_cmp = pre.cnt
    stdev_cnt_cmp = pre.std_cnt
    for filt in cnt_cmp.keys():
        if filt not in mag_cmp.keys():
            mag_cmp[filt] = np.empty(airmass.shape)
            mag_pgm[filt] = np.empty(airmass.shape)
            stdev_mag_cmp[filt] = np.empty(airmass.shape)
            stdev_mag_pgm[filt] = np.empty(airmass.shape)
        mag_cmp[filt][idx] = (cnt_cmp[filt] + post.cnt[filt]) / 2
        stdev_mag_cmp[filt][idx] = sqrt((stdev_cnt_cmp[filt]**2 + post.std_cnt[filt]**2) / 2)
        mag_pgm[filt][idx] = mid.cnt[filt]
        stdev_mag_pgm[filt][idx] = mid.std_cnt[filt]
    pre = post
    mid = None
    post = None
    idx = idx + 1

stdev_fact = 2.5 / log(10)
for filt in mag_cmp.keys():
    stdev_mag_cmp[filt] = np.abs(stdev_mag_cmp[filt] / mag_cmp[filt])
    stdev_mag_cmp[filt] = stdev_mag_cmp[filt] * stdev_fact
    mag_cmp[filt] = np.log10(mag_cmp[filt]) * -2.5
    stdev_mag_pgm[filt] = np.abs(stdev_mag_pgm[filt] / mag_pgm[filt])
    stdev_mag_pgm[filt] = stdev_mag_pgm[filt] * stdev_fact
    mag_pgm[filt] = np.log10(mag_pgm[filt]) * -2.5
    lr_cmp = linregress(airmass, mag_cmp[filt])
    print(f'CMP: k_{filt} = {lr_cmp.slope} ({lr_cmp.stderr})')
    lr_pgm = linregress(airmass, mag_pgm[filt])
    print(f'PGM: k_{filt} = {lr_pgm.slope} ({lr_pgm.stderr})')

delta_cmp = mag_cmp['B'] - mag_cmp['V']
delta_pgm = mag_pgm['B'] - mag_pgm['V']
delta_mag = delta_cmp - delta_pgm
x_delta_mag = delta_mag * airmass
delta_b = mag_cmp['B'] - mag_pgm['B']
lr_delta = linregress(x_delta_mag, delta_mag)
print(f'k"_bv = {lr_delta.slope:6.4f} {lr_delta.stderr:6.4f}')
lr_delta_b = linregress(x_delta_mag, delta_b)
print(f'k"_b = {lr_delta_b.slope:6.4f} {lr_delta_b.stderr:6.4f}')


fig, ax = plt.subplots()
#ax.plot(airmass, mag_pgm['V'], 'gx')
#ax.plot(airmass, mag_pgm['B'], 'bx')
ax.plot(x_delta_mag, x_delta_mag * lr_delta_b.slope + lr_delta_b.intercept, 'b-' )
ax.plot(x_delta_mag, delta_b, 'rx')
plt.show()