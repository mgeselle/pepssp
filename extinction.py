#!/usr/bin/python3

from argparse import ArgumentParser
from astropy.coordinates import Latitude, Longitude, EarthLocation, SkyCoord, AltAz
import astropy.units as u
from astroquery.simbad import Simbad
from dataclasses import dataclass
import numpy as np
from numpy.linalg import lstsq
from math import log10
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from typing import Sequence
from runlog import Measurement, read_measurements


@dataclass
class Star:
    pos: SkyCoord
    mag_b: float
    mag_v: float


parser = ArgumentParser(description='Determines 1st order extinction coefficients from a peppy run log.')
parser.add_argument('run_log', help='peppy run log')
parser.add_argument('--lat', help='geographic latitude of observatory',
                    default='50d6m17s')
parser.add_argument('--lon', help='geographic longitude of observatory',
                    default='8d32m12s')
parser.add_argument('--alt', type=float, help='altitude above MSL of observatory',
                    default=110)
args = parser.parse_args()

lat = Latitude(args.lat)
lon = Longitude(args.lon)
location = EarthLocation(lat=lat, lon=lon, height=args.alt)

if 'run_log' not in args:
    print('Run log not specified.')
    sys.exit(1)

run_path = Path(args.run_log)
if not run_path.is_file():
    print("Run log isn't a file or doesn't exist.")
    sys.exit(1)

measurements: Sequence[Measurement] = read_measurements(run_path)
star_names = {m.star for m in measurements}
Simbad.add_votable_fields('flux(B)', 'flux(V)', 'id(SAO)')
star_table = Simbad.query_objects(list(star_names))
stars_by_id = dict()
for s_id, ra_str, dec_str, mag_b, mag_v in star_table.iterrows('ID_SAO', 'RA', 'DEC', 'FLUX_B', 'FLUX_V'):
    pos = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
    star = Star(pos, mag_b, mag_v)
    stars_by_id[s_id] = star

airmasses = []
mags = dict()
for measurement in measurements:
    star = stars_by_id[measurement.star]
    airmasses.append(star.pos.transform_to(AltAz(obstime=measurement.time, location=location)).secz)
    for filt in measurement.cnt.keys():
        if filt not in mags:
            mags[filt] = []
        if filt == 'B':
            ref_mag = star.mag_b
        else:
            ref_mag = star.mag_v
        mags[filt].append(ref_mag + 2.5 * log10(measurement.cnt[filt]))

np_airmass = np.array(airmasses)
A = np.vstack([np_airmass, np.ones(len(np_airmass))]).T
for filt in mags.keys():
    np_mags = np.array(mags[filt])
    k, z = lstsq(A, np_mags, rcond=None)[0]
    print(f'k_{filt} = {-k:.3f}')
    plt.plot(np_airmass, np_mags, 'x')
    plt.plot(np_airmass, k * np_airmass + z, 'r')
    plt.show()
