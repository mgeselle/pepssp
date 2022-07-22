#!/usr/bin/python3
from astropy.coordinates import SkyCoord, AltAz, Latitude, Longitude, EarthLocation
import astropy.units as au
import astropy.time as atm
from astropy.utils import iers
import csv
import getopt
import math
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import lstsq
from pathlib import Path
import sys

iers.conf.auto_download = False


class CatalogStar:
    star_id = None
    pos = None
    v_mag = None
    ci_bv = None

    def __init__(self, star_id, ra, dec, pos_eq, v_mag, ci_bv):
        self.star_id = star_id
        self.pos = SkyCoord(ra, dec, unit=(au.hourangle, au.deg), equinox=atm.Time(pos_eq, format='decimalyear'))
        self.v_mag = v_mag
        self.ci_bv = ci_bv

    def airmass(self, obs_time, location):
        alt_az = self.pos.transform_to(AltAz(obstime=obs_time, location=location))
        return alt_az.secz.value


class DataPoint:
    star_id = None
    star_type = None
    time_utc = None
    airmass = None
    filter_id = None
    counts_per_sec = None

    def __init__(self, star_id, star_type, time_utc, filter_id, counts_per_sec):
        self.star_id = star_id
        self.star_type = star_type
        self.time_utc = time_utc
        self.filter_id = filter_id
        self.counts_per_sec = counts_per_sec


class Observation:
    star_id = None
    c_star_id = None
    k_star_id = None
    time_utc = None
    airmass = None
    filter_id = None
    mag_delta = None
    mag_err = None
    c_mag = None
    c_mag_err = None
    k_mag_delta = None

    def __init__(self, c_star_id, filter_id):
        self.c_star_id = c_star_id
        self.filter_id = filter_id


def load_catalog(input_file):
    result = dict()
    with open(input_file) as input_handle:
        reader = csv.DictReader(input_handle)
        for name in ('ID', 'RA', 'DEC', 'V', 'B-V', 'EQ'):
            if name not in reader.fieldnames:
                print('Missing field: {}.'.format(name))
                sys.exit(1)
        for row in reader:
            star = CatalogStar(row['ID'], row['RA'], row['DEC'], float(row['EQ']), float(row['V']), float(row['B-V']))
            result[star.star_id] = star

        return result


def load_run_log(run_log_input):
    result = dict()
    with open(run_log_input) as input_handle:
        reader = csv.DictReader(input_handle)
        for name in (
                'Timestamp', 'Index', 'StarId', 'StarType', 'IsStar', 'Filter', 'IntegrationTime', 'Count1', 'Count2',
                'Count3'):
            if name not in reader.fieldnames:
                print('Missing field: {}.'.format(name))
                sys.exit(1)
        points_by_filter = dict()
        for row in reader:
            counts = np.array((int(row['Count1']), int(row['Count2']), int(row['Count3'])))
            spread = np.max(counts) / np.min(counts)
            if spread > 1.01:
                print("Warning: line {:d}: spread of counts {:.1f}%".format(reader.line_num, (spread - 1) * 100))
            if row['IsStar'] == 'True':
                points_by_filter[row['Filter']] = DataPoint(row['StarId'],
                                                            row['StarType'],
                                                            atm.Time(row['Timestamp']),
                                                            row['Filter'],
                                                            100 * np.mean(counts) / int(row['IntegrationTime']))
            else:
                data_point = points_by_filter[row['Filter']]
                data_point.counts_per_sec = data_point.counts_per_sec - 100 * np.mean(counts) / int(
                    row['IntegrationTime'])
                if data_point.filter_id not in result:
                    result[data_point.filter_id] = []
                result[data_point.filter_id].append(data_point)
    return result


def fit_line(xvals, yvals):
    A = np.vstack([xvals, np.ones(len(xvals))]).T
    fit = lstsq(A, yvals, rcond=-1)
    a, b = fit[0]
    mse = fit[1][0] / (len(xvals) - 2)
    stderr_a = mse / math.sqrt(np.sum(np.square(xvals - np.mean(xvals))))
    stderr_b = stderr_a * math.sqrt(np.sum(np.square(xvals)))

    return a, b, stderr_a, stderr_b


def estimate_extinction(catalog, loc, obs):
    mag = []
    airmass = []
    for o in obs:
        if o.star_id not in catalog:
            print('Missing star {}'.format(o.star_id))
            continue
        airmass.append(catalog[o.star_id].airmass(o.time_utc, loc))
        std_mag = catalog[o.star_id].v_mag
        if o.filter_id == 'B':
            std_mag = std_mag + catalog[o.star_id].ci_bv
        mag.append(std_mag + 2.5 * math.log10(o.counts_per_sec))
    (k, zp, err_k, err_zp) = fit_line(np.array(airmass), np.array(mag))
    plt.plot(airmass, mag, 'r+')
    plt.plot(airmass, k * np.array(airmass) + zp, '-')
    plt.title('Primary extinction plot, {} filter'.format(obs[0].filter_id))
    plt.show()
    return k, err_k, zp


def create_observations(data_points, k_prime):
    result = []
    obs = None
    first_obs_time = None
    last_obs_time = None
    expected_types = ['CMP']
    prev_cmp = None
    curr_pgm = None
    pgm_mags = []
    cmp_mags = []
    for point in data_points:
        if point.star_type not in expected_types:
            print('Type of {} should be in {}, found {}'.format(point.star_id, ', '.join(expected_types), point.star_type))
            sys.exit(1)
        if obs is None or (point.star_type == 'CMP' and point.star_id != obs.c_star_id):
            if obs is not None:
                finish_obs(obs, pgm_mags, cmp_mags, first_obs_time, last_obs_time)
                result.append(obs)
                first_obs_time = None
                pgm_mags = []
                cmp_mags = []
            obs = Observation(point.star_id, point.filter_id)
            prev_cmp = point
            cmp_mags.append(-2.5 * math.log10(point.counts_per_sec) - k_prime * point.airmass)
            expected_types = ['PGM']
            continue
        if point.star_type == 'CMP':
            cmp_mags.append(-2.5 * math.log10(point.counts_per_sec) - k_prime * point.airmass)
            delta_m = -2.5 * math.log10(2 * curr_pgm.counts_per_sec / (prev_cmp.counts_per_sec + point.counts_per_sec))
            delta_x = curr_pgm.airmass - (prev_cmp.airmass + point.airmass) / 2
            delta_m = delta_m - k_prime * delta_x
            if curr_pgm.star_type == 'CHK':
                obs.k_mag_delta = delta_m
            else:
                pgm_mags.append(delta_m)
            prev_cmp = point
            expected_types = ['PGM', 'CHK', 'CMP']
        else:
            curr_pgm = point
            expected_types = ['CMP']
            if point.star_type == 'PGM':
                if first_obs_time is None:
                    first_obs_time = point.time_utc
                else:
                    last_obs_time = point.time_utc
                if obs.star_id is None:
                    obs.star_id = point.star_id
            elif obs.k_star_id is None:
                obs.k_star_id = point.star_id

    if obs is not None:
        finish_obs(obs, pgm_mags, cmp_mags, first_obs_time, last_obs_time)
        result.append(obs)

    return result


def finish_obs(obs, pgm_mags, cmp_mags, first_obs_time, last_obs_time):
    obs.time_utc = first_obs_time + (last_obs_time - first_obs_time) / 2
    if len(pgm_mags) > 0:
        np_mags = np.array(pgm_mags)
        obs.mag_delta = np.mean(np_mags)
        obs.mag_err = np.std(np_mags)
        np_mags = np.array(cmp_mags)
        obs.c_mag = np.mean(np_mags)
        obs.c_mag_err = np.std(np_mags)


# >>>> Main
try:
    opts, args = getopt.getopt(sys.argv[1:], 'b:c:f:o:a:h:v:')
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

catalog_file = None
run_log = None
longitude = None
latitude = None
geo_alt = 0.0
epsilon = dict()

for opt, arg in opts:
    if opt == '-c':
        catalog_file = Path(arg)
    elif opt == '-f':
        run_log = Path(arg)
    elif opt == '-o':
        longitude = Longitude(arg, unit='degree')
    elif opt == '-a':
        latitude = Latitude(arg, unit='degree')
    elif opt == '-h':
        geo_alt = float(arg)
    elif opt == '-b':
        epsilon['B'] = float(arg)
    elif opt == '-v':
        epsilon['V'] = float(arg)
    else:
        print('Unsupported option {}'.format(opt))
        sys.exit(1)

if catalog_file is None:
    print('Catalog file not specified.')
    sys.exit(1)
if not catalog_file.exists():
    print('Catalog file does not exist.')
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
print('Loading catalog from {}'.format(catalog_file.absolute()))
stars_by_id = load_catalog(catalog_file)
print('Loading run from {}'.format(run_log.absolute()))
raw_counts_by_filter = load_run_log(run_log)
for filter_counts in raw_counts_by_filter.values():
    for raw_counts in filter_counts:
        raw_counts.airmass = stars_by_id[raw_counts.star_id].airmass(raw_counts.time_utc, obs_location)

k_by_filter = {f: 0.0 for f in raw_counts_by_filter}
obs_by_filter = dict()
for filt_id in k_by_filter:
    extinction_pts = [p for p in raw_counts_by_filter[filt_id] if p.star_type == 'EXT']
    if len(extinction_pts) > 0:
        print('Estimating primary extinction for filter {}'.format(filt_id))
        (k, err_k, zp) = estimate_extinction(stars_by_id, obs_location, extinction_pts)
        print("k'({}) = {:5.3f}({:5.3f}), zeta = {:6.3f}".format(filt_id, -k, err_k, zp))
        k_by_filter[filt_id] = -k
    else:
        k_by_filter[filt_id] = 0.0

    obs_pts = [p for p in raw_counts_by_filter[filt_id] if p.star_type in ['CMP', 'CHK', 'PGM']]
    obs_by_filter[filt_id] = create_observations(obs_pts, k_by_filter[filt_id])
    if filt_id not in epsilon:
        epsilon[filt_id] = 0.0
    print('Using epsilon_{} = {:5.3f}'.format(filt_id, epsilon[filt_id]))
    for obs in obs_by_filter[filt_id]:
        if filt_id == 'V':
            mag = stars_by_id[obs.c_star_id].v_mag + obs.mag_delta
        else:
            mag = stars_by_id[obs.c_star_id].ci_bv + stars_by_id[obs.c_star_id].v_mag + obs.mag_delta

        mag = mag + epsilon[filt_id] * stars_by_id[obs.star_id].ci_bv
        c_mag = obs.c_mag + zp + epsilon[filt_id] * stars_by_id[obs.c_star_id].ci_bv

        print('Observation {}@{}'.format(obs.star_id, filt_id))
        print('Mag: {:5.3f}({:5.3f})'.format(mag, obs.mag_err))
        print('CMag: {:5.3f}({:5.3f}), Abs: {:5.3f}'.format(obs.c_mag, obs.c_mag_err, c_mag))
        print('KMag: {:5.3f}'.format(obs.k_mag_delta + obs.c_mag))
        print('Airmass: {:5.3f}'.format(stars_by_id[obs.star_id].airmass(obs.time_utc, obs_location)))
        obs_time_hc = obs.time_utc.utc + obs.time_utc.light_travel_time(stars_by_id[obs.star_id].pos,
                                                                     kind='heliocentric',
                                                                     location=obs_location)
        print('HJD: {:14.6f}'.format(obs_time_hc.jd))
