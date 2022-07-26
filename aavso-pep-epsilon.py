#!/usr/bin/python3

from astropy.coordinates import SkyCoord, Longitude, Latitude, EarthLocation, AltAz
import astropy.time as atime
from dataclasses import dataclass
import getopt
from math import sqrt, log, fabs
import numpy as np
import sys
from typing import Union, Sequence

from runlog import Measurement, read_measurements


@dataclass
class CalibPair:
    loc_red: SkyCoord
    loc_blue: SkyCoord
    delta_v: float
    delta_b: float


pairs_by_cmp = {
    'SAO 65890': CalibPair(SkyCoord('17h15m02.8343634s', '+36d48m32.98434s', frame='icrs'),
                           SkyCoord('17h17m40.253198052s', '+37d17m29.422635288s', frame='icrs'),
                           1.465, 0.071)
}


def calc_airmass(pos: SkyCoord, time: atime.Time, loc: EarthLocation) -> float:
    return pos.transform_to(AltAz(obstime=time, location=loc)).secz


def evaluate_measurements(measurements: Sequence[Measurement], location: EarthLocation, k_b: Union[None, float]):
    if measurements[0].star not in pairs_by_cmp:
        raise ValueError(f'No pair for comparison {measurements[0].star}')
    pair = pairs_by_cmp[measurements[0].star]

    pre = None
    post = None
    blue = None

    num_blue = int((len(measurements) - 1) / 2)
    airmass = np.empty(num_blue)
    mag_vis_red = np.empty(num_blue)
    mag_vis_blue = np.empty(num_blue)
    stdev_mag_vis_red = np.empty(num_blue)
    stdev_mag_vis_blue = np.empty(num_blue)
    if 'B' in measurements[0].cnt and k_b is not None:
        mag_blue_red = np.empty(num_blue)
        mag_blue_blue = np.empty(num_blue)
        stdev_mag_blue_red = np.empty(num_blue)
        stdev_mag_blue_blue = np.empty(num_blue)
    else:
        mag_blue_red = None
        mag_blue_blue = None
        stdev_mag_blue_red = None
        stdev_mag_blue_blue = None
    idx = 0
    for measurement in measurements:
        if pre is None:
            pre = measurement
            continue
        elif blue is None:
            blue = measurement
            continue
        elif post is None:
            post = measurement

        pre_airmass = calc_airmass(pair.loc_red, pre.time, location)
        blue_airmass = calc_airmass(pair.loc_blue, blue.time, location)
        post_airmass = calc_airmass(pair.loc_red, post.time, location)
        airmass[idx] = (pre_airmass + blue_airmass + post_airmass) / 3

        mag_vis_red[idx] = (pre.cnt['V'] + post.cnt['V']) / 2
        stdev_mag_vis_red[idx] = sqrt((pre.std_cnt['V']**2 + post.std_cnt['V']**2) / 4)
        mag_vis_blue[idx] = blue.cnt['V']
        stdev_mag_vis_blue[idx] = blue.std_cnt['V']
        if mag_blue_red is not None:
            mag_blue_red[idx] = (pre.cnt['B'] + post.cnt['B']) / 2
            stdev_mag_blue_red[idx] = sqrt((pre.std_cnt['B']**2 + post.std_cnt['B']**2) / 4)
            mag_blue_blue[idx] = blue.cnt['B']
            stdev_mag_blue_blue[idx] = blue.std_cnt['B']

        pre = post
        blue = None
        post = None
        idx = idx + 1

    delta_std_v = pair.delta_v
    delta_std_bv = pair.delta_b - pair.delta_v

    factor = 2.5 / log(10)
    stdev_mag_vis_red = stdev_mag_vis_red / mag_vis_red * factor
    mag_vis_red = np.log10(mag_vis_red) * -2.5
    stdev_mag_vis_blue = stdev_mag_vis_blue / mag_vis_blue * factor
    mag_vis_blue = np.log10(mag_vis_blue) * -2.5

    delta_v = np.mean(mag_vis_blue - mag_vis_red)
    stdev_delta_v = sqrt((np.mean(stdev_mag_vis_blue**2 + stdev_mag_vis_red**2)) / num_blue)

    epsilon_v = (delta_std_v - delta_v) / delta_std_bv
    stdev_epsilon_v = fabs(stdev_delta_v / delta_std_bv)
    print(f'epsilon_v = {epsilon_v:6.4f} ({stdev_epsilon_v:6.4f})')

    if mag_blue_red is not None:
        stdev_mag_blue_red = stdev_mag_blue_red / mag_blue_red * factor
        mag_blue_red = np.log10(mag_blue_red) * -2.5
        stdev_mag_blue_blue = stdev_mag_blue_blue / mag_blue_blue * factor
        mag_blue_blue = np.log10(mag_blue_blue) * -2.5

        delta_b = np.mean(mag_blue_blue - mag_blue_red - airmass * delta_std_bv * k_b)
        stdev_delta_b = sqrt((np.mean(stdev_mag_blue_blue**2) + np.mean(stdev_mag_blue_red**2)) / num_blue)

        epsilon_b = (pair.delta_b - delta_b) / delta_std_bv
        stdev_epsilon_b = fabs(stdev_delta_b / delta_std_bv)
        print(f'epsilon_b = {epsilon_b:6.4f} ({stdev_epsilon_b:6.4f})')


def usage_exit():
    print(f'Usage: {sys.argv[0]} -o longitude -a latitude [-h altitude] [-s k"_b] -f PEP_run_log')
    sys.exit(1)


try:
    opts, args = getopt.getopt(sys.argv[1:], 'f:o:a:h:s:')
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

run_log = None
longitude = None
latitude = None
geo_alt = 0.0
k_b = None

for opt, arg in opts:
    if opt == '-f':
        run_log = arg
    elif opt == '-o':
        longitude = Longitude(arg)
    elif opt == '-a':
        latitude = Latitude(arg)
    elif opt == '-h':
        geo_alt = float(arg)
    elif opt == '-s':
        k_b = float(arg)
    else:
        print('Unsupported option {}'.format(opt))
        sys.exit(1)

if run_log is None:
    print('Run log not specified.')
    sys.exit(1)
if longitude is None:
    print('Longitude not specified.')
    sys.exit(1)
if latitude is None:
    print('Latitude not specified')
    sys.exit(1)

obs_location = EarthLocation(longitude, latitude, geo_alt)

loaded = read_measurements(run_log)
evaluate_measurements(loaded, obs_location, k_b)