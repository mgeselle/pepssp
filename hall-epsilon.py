#!/usr/bin/python3
"""Determines the visual transformation coefficient epsilon_v by the method detailed in
Hall, D.D Blue-Red Star Pairs for Determining Transformation Coefficients"""
from astropy.coordinates import SkyCoord, AltAz, Latitude, Longitude, EarthLocation
import astropy.units as au
import astropy.time as atm
from astropy.utils import iers
import csv
import getopt
import math
import sys
from pathlib import Path
import numpy as np


iers.conf.auto_download = False

# These are the AAVSO red-pairs. idx is the pair ID, role is the role of the
# star in the pair. The star names given in https://www.aavso.org/obtaining-your-pep-epsilonv-coefficient
# have been replaced by SAO catalog IDs. This is because the Gemini 1 in my G11 mount has the SAO
# catalog in its built-in object database. The pair ID is an index into the pairs array below.
pairs_by_star = {
    'SAO 37434': {'idx': 0, 'role': 'r'},
    'SAO 37418': {'idx': 0, 'role': 'b'},
    'SAO 38890': {'idx': 1, 'role': 'r'},
    'SAO 38893': {'idx': 1, 'role': 'b'},
    'SAO 112098': {'idx': 2, 'role': 'r'},
    'SAO 112096': {'idx': 2, 'role': 'b'},
    'SAO 62019': {'idx': 3, 'role': 'r'},
    'SAO 62010': {'idx': 3, 'role': 'b'},
    'SAO 121157': {'idx': 4, 'role': 'r'},
    'SAO 121170': {'idx': 4, 'role': 'b'},
    'SAO 65890': {'idx': 5, 'role': 'r'},
    'SAO 65921': {'idx': 5, 'role': 'b'},
    'SAO 122686': {'idx': 6, 'role': 'r'},
    'SAO 122687': {'idx': 6, 'role': 'b'}
}

# Red-blue pairs by pair ID.
# cstl: constellation
# rvmag: red star's V magnitude
# bvmag: blue star's V magnitude
# delta_bv: difference in (B - V), blue star minus red star
# r_pos: position of red star
# b_pos: position of blue star
pairs = [
    {'cstl': 'And', 'rvmag': 4.955, 'bvmag': 4.949, 'delta_bv': -0.709,
     'r_pos': SkyCoord('01:42:00', '42:37:00', unit=(au.hourangle, au.deg)),
     'b_pos': SkyCoord('01:41:00', '40:35:00', unit=(au.hourangle, au.deg))},
    {'cstl': 'Per', 'rvmag': 4.357, 'bvmag': 5.833, 'delta_bv': -1.406,
     'r_pos': SkyCoord('03:31:00', '48:00:00', unit=(au.hourangle, au.deg)),
     'b_pos': SkyCoord('03:31:00', '48:06:00', unit=(au.hourangle, au.deg))},
    {'cstl': 'Ori', 'rvmag': 6.03, 'bvmag': 7.32, 'delta_bv': -1.245,
     'r_pos': SkyCoord('04:49:00', '03:35:00', unit=(au.hourangle, au.deg)),
     'b_pos': SkyCoord('04:49:00', '03:59:00', unit=(au.hourangle, au.deg))},
    {'cstl': 'LMi', 'rvmag': 5.50, 'bvmag': 5.878, 'delta_bv': -1.03,
     'r_pos': SkyCoord('10:24:00', '33:43:00', unit=(au.hourangle, au.deg)),
     'b_pos': SkyCoord('10:23:00', '33:54:00', unit=(au.hourangle, au.deg))},
    {'cstl': 'Ser', 'rvmag': 2.65, 'bvmag': 5.58, 'delta_bv': -1.13,
     'r_pos': SkyCoord('15:44:00', '06:26:00', unit=(au.hourangle, au.deg)),
     'b_pos': SkyCoord('15:45:00', '05:27:00', unit=(au.hourangle, au.deg))},
    {'cstl': 'Her', 'rvmag': 3.161, 'bvmag': 4.626, 'delta_bv': -1.394,
     'r_pos': SkyCoord('17:15:00', '36:49:00', unit=(au.hourangle, au.deg)),
     'b_pos': SkyCoord('17:18:00', '37:18:00', unit=(au.hourangle, au.deg))},
    {'cstl': 'Oph', 'rvmag': 7.805, 'bvmag': 8.315, 'delta_bv': -1.23,
     'r_pos': SkyCoord('17:44:00', '05:15:00', unit=(au.hourangle, au.deg)),
     'b_poa': SkyCoord('17:44:00', '05:43:00', unit=(au.hourangle, au.deg))}
]


def usage_exit():
    print('Usage: {} -o longitude -a latitude [-s k"_b] -f PEP_run_log'.format(sys.argv[0]))
    sys.exit(1)


def finish_item(star_idx, deflection, red_measurements, blue_measurements):
    if star_idx['role'] == 'b':
        blue_measurements.append(deflection)
    else:
        red_measurements.append(deflection)
    print(star_idx['role'] + '-> ', end=' ')
    for filt in sorted(deflection.keys()):
        print('{}={:4.0f}'.format(filt, deflection[filt]), end=' ')
    print()


def load_run(run_log):
    """Loads a PEP run from a csv file containing one line per measurement. A measurement
    is either for a star of for the sky background and always has 3 deflections. The file is
    expected to have column headers in the first line which determine which field is where.
    Expected fields:

    Index: index of the measurement. A pair of star-sky measurements shares a common index.
    StarId: ID of the star, one of the SAO catalog IDs from the pairs_by_star dictionary.
    IsStar: True or False: flag which indicates whether a row is a star or sky background.
    Filter: Filter used for measurement (e.g.: V).
    IntegrationTime: integration time in hundredths of a second.
    Count1, Count2, Count3: the individual deflections for this measurement as read from the photometer.

    Note: in theory this program could also handle other/multiple filters. This hasn't been implemented
    (yet), because use of other filters requires the second order extinction coefficient to be
    known."""
    red_measurements = []
    blue_measurements = []
    l_pair_id = None
    item_idx = -1
    item = {}
    with open(run_log) as in_handle:
        reader = csv.DictReader(in_handle)
        for name in ('Index', 'StarId', 'StarType', 'IsStar', 'Filter', 'IntegrationTime',
                     'Count1', 'Count2', 'Count3'):
            if name not in reader.fieldnames:
                print('Missing field: {}'.format(name))
                sys.exit(1)
        star_index = None
        start_time = None
        end_time = None
        for row in reader:
            if row['StarId'] not in pairs_by_star:
                continue
            if item_idx != int(row['Index']):
                if item:
                    finish_item(star_index, item, red_measurements, blue_measurements)
                    item = {}
                item_idx = int(row['Index'])
            star_index = pairs_by_star[row['StarId']]
            if not l_pair_id and l_pair_id != 0:
                l_pair_id = star_index['idx']
                start_time = atm.Time(row['Timestamp'])
                print('Computing epsilon_v from pair in {}.'.format(pairs[l_pair_id]['cstl']))
                print('Deflections in counts/s:')
            elif star_index['idx'] != l_pair_id:
                continue
            end_time = atm.Time(row['Timestamp'])
            counts_per_sec = (int(row['Count1']) + int(row['Count2']) + int(row['Count3'])) / 3
            counts_per_sec = counts_per_sec / int(row['IntegrationTime'])
            counts_per_sec = counts_per_sec * 100
            if row['IsStar'] == 'True':
                item[row['Filter']] = counts_per_sec
            else:
                item[row['Filter']] = item[row['Filter']] - counts_per_sec
        if item:
            finish_item(star_index, item, red_measurements, blue_measurements)
    return l_pair_id, start_time, end_time, red_measurements, blue_measurements


def calc_airmass(star_pos, loc, obs_tm):
    alt_az = star_pos.transform_to(AltAz(obstime=obs_tm, location=loc))
    return alt_az.secz.value


try:
    opts, args = getopt.getopt(sys.argv[1:], 'f:o:a:h:s:')
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

catalog_file = None
run_log = None
longitude = None
latitude = None
geo_alt = 0.0
# Assumed value for secondary bv extinction coefficient
k_bv = -0.04

for opt, arg in opts:
    if opt == '-f':
        run_log = Path(arg)
    elif opt == '-o':
        longitude = Longitude(arg, unit='degree')
    elif opt == '-a':
        latitude = Latitude(arg, unit='degree')
    elif opt == '-h':
        geo_alt = float(arg)
    elif opt == '-s':
        k_bv = float(arg)
    else:
        print('Unsupported option {}'.format(opt))
        sys.exit(1)

if run_log is None:
    print('Run log not specified.')
    sys.exit(1)
if not run_log.exists():
    print('Run log does not exist.')
    sys.exit(1)
if longitude is None:
    print('Longitude not specified.')
    sys.exit(1)
if latitude is None:
    print('Latitude not specified')
    sys.exit(1)

obs_location = EarthLocation(longitude, latitude, geo_alt)

pair_id, start_tm, end_tm, red, blue = load_run(run_log)
mags = {}
for k in red[0]:
    mags[k] = []
for idx, blue_defl in enumerate(blue):
    for k in blue_defl:
        mag = 2 * blue_defl[k] / (red[idx][k] + red[idx + 1][k])
        mags[k].append(-2.5 * math.log10(mag))

pair = pairs[pair_id]
epsilon_by_filter = dict()
for k in mags:
    np_mags = np.array(mags[k])
    delta_m0 = np.mean(np_mags)
    sigma_m0 = np.std(np_mags)
    print('Filter: {}, delta_m0 = {:7.4f}({:7.4f})'.format(k, delta_m0, sigma_m0))
    if k == 'V':
        epsilon = (pair['bvmag'] - pair['rvmag'] - delta_m0) / pair['delta_bv']
        print('Epsilon_v = {:7.4f}'.format(epsilon))
    else:
        airmass = calc_airmass(pair['r_pos'], obs_location, start_tm)
        airmass = airmass + calc_airmass(pair['b_pos'], obs_location, start_tm)
        airmass = airmass + calc_airmass(pair['r_pos'], obs_location, end_tm)
        airmass = airmass + calc_airmass(pair['b_pos'], obs_location, end_tm)
        airmass = airmass / 4
        delta_b = pair['delta_bv'] + pair['bvmag'] - pair['rvmag']
        epsilon = (delta_b - delta_m0 + k_bv * airmass * pair['delta_bv']) / pair['delta_bv']
        print('Epsilon_b = {:7.4f}'.format(epsilon))
    epsilon_by_filter[k] = epsilon

if 'B' in epsilon_by_filter and 'V' in epsilon_by_filter:
    print('Mu_bv = {:7.4f}'.format(1 / (1 - epsilon_by_filter['B'] + epsilon_by_filter['V'])))


