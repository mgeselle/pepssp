#!/usr/bin/python3
import sys
from argparse import ArgumentParser
import astropy.time as atime
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, Latitude, Longitude
from dataclasses import dataclass
from math import log, log10, sqrt
import numpy as np
from os import PathLike
from pathlib import Path
from typing import Union, Dict, Literal, Sequence, SupportsFloat

import aavsosp
from runlog import Measurement, read_measurements


@dataclass
class Magnitude:
    mag: Union[float | SupportsFloat]
    mag_err: float


@dataclass()
class MeasuredStar:
    name: str
    mags: Dict[Literal['B', 'V', 'R', 'I'], Magnitude]


@dataclass
class MeasuredEnsemble:
    time: atime.Time
    pgm: MeasuredStar
    cmp: MeasuredStar
    chk: MeasuredStar


@dataclass
class TransformedEnsemble:
    time: atime.Time
    pgm: MeasuredStar
    chk: MeasuredStar
    chk_i: MeasuredStar
    cmp_i: MeasuredStar
    chk_ref: Dict[Literal['B', 'V', 'R', 'I'], float]
    pgm_x: float
    cmp_x: float
    chk_x: float


@dataclass
class TransformParams:
    starparm: Union[str, bytes, PathLike]
    epsilon: Dict[Literal['B', 'V'], float]
    k: Dict[Literal['B', 'V'], float]
    kk: float
    location: EarthLocation
    obs_code: str
    scope: str


def form_ensembles(measurements: Sequence[Measurement]) -> Sequence[MeasuredEnsemble]:
    result = []
    cmp_counts = dict()
    cmp_count_errs = dict()
    last_cmp_count = dict()
    last_cmp_count_err = dict()
    pgm_counts = dict()
    pgm_count_errs = dict()
    cmp_idx = 0
    pgm_name = None
    cmp_name = None
    chk_name = None
    chk_count = dict()
    chk_count_err = dict()
    time = None
    err_factor = 2.5 / log(10)
    for measurement in measurements:
        if measurement.star_type == 'CMP':
            if cmp_name is None:
                cmp_name = measurement.star
            for filt in measurement.cnt.keys():
                if filt not in cmp_counts:
                    cmp_counts[filt] = np.empty(5)
                    cmp_count_errs[filt] = np.empty(5)
                cmp_count_errs[filt][cmp_idx] = measurement.std_cnt[filt]
                cmp_counts[filt][cmp_idx] = measurement.cnt[filt]
                last_cmp_count[filt] = measurement.cnt[filt]
                last_cmp_count_err[filt] = measurement.std_cnt[filt]
            cmp_idx = cmp_idx + 1
        elif measurement.star_type == 'PGM':
            if pgm_name is None:
                pgm_name = measurement.star
            if cmp_idx == 0:
                # Back-to-back measurement of same star...
                for filt in last_cmp_count.keys():
                    cmp_counts[filt][0] = last_cmp_count[filt]
                    cmp_count_errs[filt][0] = last_cmp_count_err[filt]
                cmp_idx = 1
            if cmp_idx == 2:
                time = measurement.time
            for filt in measurement.cnt.keys():
                if filt not in pgm_counts:
                    pgm_counts[filt] = np.empty(3)
                    pgm_count_errs[filt] = np.empty(3)
                pgm_count_errs[filt][cmp_idx - 1] = measurement.std_cnt[filt]
                pgm_counts[filt][cmp_idx - 1] = measurement.cnt[filt]
        elif measurement.star_type == 'CHK':
            if chk_name is None:
                chk_name = measurement.star
            for filt in measurement.cnt.keys():
                chk_count[filt] = measurement.cnt[filt]
                chk_count_err[filt] = measurement.std_cnt[filt]

        if cmp_idx == 5:
            cmp_mags = dict()
            pgm_mags = dict()
            chk_mags = dict()
            for filt in cmp_counts.keys():
                cmp_count_errs[filt] = cmp_count_errs[filt] / cmp_counts[filt] * err_factor
                cmp_counts[filt] = np.log10(cmp_counts[filt]) * -2.5
                cmp_mags[filt] = Magnitude(np.mean(cmp_counts[filt]), sqrt(np.mean(cmp_count_errs[filt] ** 2) / 4))
                pgm_count_errs[filt] = pgm_count_errs[filt] / pgm_counts[filt] * err_factor
                pgm_counts[filt] = np.log10(pgm_counts[filt]) * -2.5
                pgm_mags[filt] = Magnitude(np.mean(pgm_counts[filt]), sqrt(np.mean(pgm_count_errs[filt] ** 2) / 3))
                chk_mags[filt] = Magnitude(-2.5 * log10(chk_count[filt]),
                                           err_factor * chk_count_err[filt] / chk_count[filt])
            cmp_star = MeasuredStar(cmp_name, cmp_mags)
            pgm_star = MeasuredStar(pgm_name, pgm_mags)
            chk_star = MeasuredStar(chk_name, chk_mags)
            ensemble = MeasuredEnsemble(time, pgm_star, cmp_star, chk_star)
            result.append(ensemble)

            cmp_idx = 0
            cmp_name = None
            pgm_name = None
            chk_name = None
            time = None
            # Clearing out filters, because they might change from
            # observation to observation.
            cmp_counts = dict()
            cmp_count_errs = dict()
            pgm_counts = dict()
            pgm_count_errs = dict()
            chk_count = dict()
            chk_count_err = dict()

    return result


def _calc_airmass(pos: SkyCoord, loc: EarthLocation, time: atime.Time):
    return pos.transform_to(AltAz(obstime=time, location=loc)).secz


def transform(observations: Sequence[MeasuredEnsemble], params: TransformParams) -> Sequence[TransformedEnsemble]:
    result = []
    for observation in observations:
        ensemble = aavsosp.locate_ensemble(params.starparm, observation.pgm.name)
        if ensemble is None:
            print(f'Star {observation.pgm.name} not found in starparm file.')
            continue
        mags = dict()
        chk_mags = dict()
        chk_i_mags = dict()
        m_chk_ref = dict()
        pgm_airmass = _calc_airmass(ensemble.pgm.pos, params.location, observation.time)
        cmp_airmass = _calc_airmass(ensemble.cmp.pos, params.location, observation.time)
        chk_airmass = _calc_airmass(ensemble.chk.pos, params.location, observation.time)
        for filt in observation.pgm.mags.keys():
            if filt == 'B':
                m_ref = ensemble.cmp.bmv + ensemble.cmp.vmag
                if ensemble.chk.bmv is None:
                    m_chk_ref[filt] = None
                else:
                    m_chk_ref[filt] = ensemble.chk.bmv + ensemble.chk.vmag
            elif filt == 'V':
                m_ref = ensemble.cmp.vmag
                m_chk_ref[filt] = ensemble.chk.vmag
            else:
                print(f'Unsupported filter: {filt}')
                continue
            pgm_mag = m_ref + observation.pgm.mags[filt].mag - observation.cmp.mags[filt].mag + \
                params.epsilon[filt] * (ensemble.pgm.bmv - ensemble.cmp.bmv) - \
                params.k[filt] * (pgm_airmass - cmp_airmass)
            chk_mag = m_ref + observation.chk.mags[filt].mag - observation.cmp.mags[filt].mag - \
                params.k[filt] * (chk_airmass - cmp_airmass)
            chk_i_mag = chk_mag
            if ensemble.chk.bmv is not None:
                chk_mag = chk_mag + params.epsilon[filt] * (ensemble.chk.bmv - ensemble.cmp.bmv)
            chk_mag = chk_mag + params.k[filt] * (chk_airmass - cmp_airmass)
            if filt == 'B':
                pgm_cmp_am = (pgm_airmass + cmp_airmass) / 2
                chk_cmp_am = (chk_airmass + cmp_airmass) / 2
                bmv_pgm = observation.pgm.mags['B'].mag - observation.pgm.mags['V'].mag
                bmv_cmp = observation.cmp.mags['B'].mag - observation.cmp.mags['V'].mag
                bmv_chk = observation.chk.mags['B'].mag - observation.chk.mags['V'].mag
                pgm_mag = pgm_mag - params.kk * pgm_cmp_am * (bmv_pgm - bmv_cmp)
                chk_mag = chk_mag - params.kk * chk_cmp_am * (bmv_chk - bmv_cmp)
                chk_i_mag = chk_i_mag - params.kk * chk_cmp_am * (bmv_chk - bmv_cmp)
                err_pgm = sqrt(((1 - params.kk * pgm_cmp_am) * observation.pgm.mags['B'].mag_err) ** 2 +
                               ((1 - params.kk * pgm_cmp_am) * observation.cmp.mags['B'].mag_err) ** 2 +
                               (params.kk * pgm_cmp_am * observation.pgm.mags['V'].mag_err) ** 2 +
                               (params.kk * pgm_cmp_am * observation.cmp.mags['V'].mag_err) ** 2)
                err_chk = sqrt(((1 - params.kk * chk_cmp_am) * observation.chk.mags['B'].mag_err) ** 2 +
                               ((1 - params.kk * chk_cmp_am) * observation.cmp.mags['B'].mag_err) ** 2 +
                               (params.kk * chk_cmp_am * observation.chk.mags['V'].mag_err) ** 2 +
                               (params.kk * chk_cmp_am * observation.cmp.mags['V'].mag_err) ** 2)
            else:
                err_pgm = sqrt(observation.pgm.mags[filt].mag_err ** 2 + observation.cmp.mags[filt].mag_err ** 2)
                err_chk = sqrt(observation.chk.mags[filt].mag_err ** 2 + observation.cmp.mags[filt].mag_err ** 2)
            mags[filt] = Magnitude(pgm_mag, err_pgm)
            chk_mags[filt] = Magnitude(chk_mag, err_chk)
            chk_i_mags[filt] = Magnitude(chk_i_mag, err_chk)
        pgm = MeasuredStar(ensemble.pgm.usename, mags)
        chk = MeasuredStar(ensemble.chk.usename, chk_mags)
        chk_i = MeasuredStar(ensemble.chk.usename, chk_i_mags)
        cmp_i = MeasuredStar(ensemble.cmp.usename, observation.cmp.mags)
        xf = TransformedEnsemble(observation.time, pgm, chk, chk_i, cmp_i, m_chk_ref,
                                 pgm_airmass, cmp_airmass, chk_airmass)
        result.append(xf)

    return result


def report(stars: Sequence[TransformedEnsemble]):
    if len(stars) == 0:
        return
    print()
    print('++++++++++++++++ Results +++++++++++++')
    print()
    for star in stars:
        usename = star.pgm.name
        time_str = star.time.strftime('%d.%m.%Y %H:%M:%S')
        print(f'{usename} at {time_str} UTC (JD={star.time.jd:13.5f})')
        print()
        print('  PGM   ERR   CHK   ERR   CHK Ref')
        print('---------------------------------')
        for filt in star.pgm.mags:
            m_pgm = star.pgm.mags[filt].mag
            me_pgm = star.pgm.mags[filt].mag_err
            m_chk = star.chk.mags[filt].mag
            me_chk = star.chk.mags[filt].mag_err
            mr_chk = star.chk_ref[filt]
            if mr_chk is not None:
                print(f'{filt} {m_pgm:5.3f} {me_pgm:5.3f} {m_chk:5.3f} {me_chk:5.3f} {mr_chk:5.3f}')
            else:
                print(f'{filt} {m_pgm:5.3f} {me_pgm:5.3f} {m_chk:5.3f} {me_chk:5.3f} -----')


def write_aavso(stars: Sequence[TransformedEnsemble], file_name: Union[str, bytes, PathLike],
                params: TransformParams):
    if len(stars) == 0:
        return
    file_path = Path(file_name)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    with file_path.open('w') as out:
        print('#TYPE=Extended', file=out)
        print(f'#OBSCODE={params.obs_code}', file=out)
        print('#SOFTWARE=pep-reduce.py, see: https://github.com/mgeselle/pepssp', file=out)
        print('#DELIM=,', file=out)
        print('#DATE=JD', file=out)
        print('#OBSTYPE=PEP', file=out)

        lat_i = int(round(params.location.lat.deg))
        if lat_i > 0:
            lat_ns = 'N'
        else:
            lat_ns = 'S'
        lon_i = int(round(params.location.lon.deg))
        if lon_i > 0:
            lon_ew = 'E'
        else:
            lon_ew = 'W'
        notes_pfx = f'SCOPE={params.scope}|SENSOR=SSP3|LOC={lat_i}{lat_ns}/{lon_i}{lon_ew}|INDEX=BV'

        for star in stars:
            usename = star.pgm.name
            pgm_x = star.pgm_x
            cmp_x = star.cmp_x
            chk_x = star.chk_x
            jd = star.time.jd
            for filt in star.pgm.mags.keys():
                notes = notes_pfx + f'|CX={cmp_x:5.3f}|KX={chk_x:5.3f}|K_{filt}={params.k[filt]:4.2f}'
                if filt == 'B':
                    notes = notes + f'|KK_B={params.kk:5.2f}|TB_BV={params.epsilon[filt]:4.2f}'
                elif filt == 'V':
                    notes = notes + f'|TV_BV={params.epsilon[filt]:4.2f}'
                pgm_mag = star.pgm.mags[filt].mag
                pgm_err = star.pgm.mags[filt].mag_err
                cname = star.cmp_i.name
                cmag = star.cmp_i.mags[filt].mag
                kname = star.chk_i.name
                kmag = star.chk_i.mags[filt].mag
                line = f'{usename},{jd:13.5f},{pgm_mag:5.3f},{pgm_err:5.3f},{filt}'
                line = line + f',YES,STD,{cname},{cmag:5.3f},{kname},{kmag:5.3f}'
                line = line + f',{pgm_x:5.3f},na,PEP,{notes}'

                print(line, file=out)


parser = ArgumentParser(description='Perform data reduction on a peppy run log.')
parser.add_argument('run_log', help='peppy run log')
parser.add_argument('-o', '--out', help='output AAVSO photometry file')
parser.add_argument('--kb', type=float, help='primary extinction coefficient for B',
                    default=0.351)
parser.add_argument('--kv', type=float, help='primary extinction coefficient for V',
                    default=0.201)
parser.add_argument('--kk', type=float, help='secondary extinction coefficient for B',
                    default=-0.026)
parser.add_argument('--eb', type=float, help='transformation coefficient for B',
                    default=0.1476)
parser.add_argument('--ev', type=float, help='transformation coefficient for V',
                    default=-0.0196)
parser.add_argument('--starparm', help='AAVSO starparm file',
                    default=str(Path.home() / 'Documents/astro/starparm_2021Aug27.csv'))
parser.add_argument('--lat', help='geographic latitude of observatory',
                    default='50d6m17s')
parser.add_argument('--lon', help='geographic longitude of observatory',
                    default='8d32m12s')
parser.add_argument('--alt', type=float, help='altitude above MSL of observatory',
                    default=110)
parser.add_argument('--obs', help='AAVSO Observer Code', default='GMV')
parser.add_argument('--sco', help='Telescope used', default='8IN_RC')

args = parser.parse_args()

eps = {'B': args.eb, 'V': args.ev}
k = {'B': args.kb, 'V': args.kv}
lat = Latitude(args.lat)
lon = Longitude(args.lon)

# noinspection PyTypeChecker
xform_params = TransformParams(args.starparm, eps, k, args.kk,
                               EarthLocation(lat=lat, lon=lon, height=args.alt),
                               args.obs, args.sco)
if 'run_log' not in args:
    print('Run log not specified.')
    sys.exit(1)

m = read_measurements(args.run_log)
e = form_ensembles(m)
x = transform(e, xform_params)
report(x)
if 'out' in args and args.out is not None:
    write_aavso(x, args.out, xform_params)
