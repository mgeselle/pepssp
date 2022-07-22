#!/usr/bin/python3

from astropy.coordinates import EarthLocation, Latitude, Longitude, AltAz, SkyCoord
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from astroquery.simbad import Simbad
from dataclasses import dataclass
import getopt
import json
from math import fabs
from os import PathLike
from pathlib import Path
import sys
from typing import Sequence, Tuple, Union

import peprun


@dataclass
class ExtStar:
    id: str
    sao: str
    mag_b: float
    mag_v: float
    rising: bool
    airmass: float


def get_ra_limits(time: Time, location: EarthLocation, dec_limit: float) -> Tuple[Longitude, Longitude]:
    siderial_t = time.sidereal_time('apparent', longitude=location.lon).deg
    airmass = 1
    ra_angle = siderial_t
    while airmass < 3:
        ra_angle = ra_angle - 5.0
        if ra_angle < 0:
            ra_angle = ra_angle + 360.0
        pos = SkyCoord(ra_angle * u.deg, dec_limit * u.deg)
        airmass = calc_airmass(pos, time, location)
    setting = Longitude(ra_angle * u.deg)
    rising_angle = siderial_t + (siderial_t - ra_angle)
    if rising_angle > 360.0:
        rising_angle = rising_angle - 360.0
    rising = Longitude(rising_angle * u.deg)

    return setting, rising


def calc_airmass(pos: SkyCoord, time: Time, location: EarthLocation) -> float:
    return pos.transform_to(AltAz(obstime=time, location=location)).secz


def query_simbad_objects(setting: Longitude, rising: Longitude, dec_limit: float) -> Table:
    query = f"dec > {dec_limit} & spstring ~ '^A' & Vmag >= 3 & Vmag < 6"
    query = query + f" & (maintype = '*' | maintype = 'PM*')"
    setting_ha = setting.deg
    rising_ha = rising.deg
    if rising_ha > setting_ha:
        query = query + f" & ra > {setting_ha:5.1f} & ra < {rising_ha:5.1f}"
    else:
        query = query + f" & (ra > {setting_ha:5.1f} | ra < {rising_ha:5.1f})"

    print(f'Submitting query: "{query}"')
    Simbad.add_votable_fields('flux(B)', 'flux(V)', 'id(SAO)')
    return Simbad.query_criteria(query)


def pick_stars(candidates: Table, time: Time, location: EarthLocation) -> Sequence[ExtStar]:
    all_stars_rising = []
    all_stars_setting = []
    # noinspection PyTypeChecker
    for s_id, ra_str, dec_str, mag_b, mag_v, sao in \
            candidates.iterrows('MAIN_ID', 'RA', 'DEC', 'FLUX_B', 'FLUX_V', 'ID_SAO'):
        if fabs(mag_b - mag_v) > 0.2 or sao == '':
            continue
        pos = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
        altaz = pos.transform_to(AltAz(obstime=time, location=location))
        star = ExtStar(s_id, sao, mag_b, mag_v, altaz.az.deg < 180, altaz.secz.value)
        if star.rising:
            all_stars_rising.append(star)
        else:
            all_stars_setting.append(star)

    all_stars_rising.sort(key=lambda x: -x.airmass)
    all_stars_setting.sort(key=lambda x: -x.airmass)
    range_rising = all_stars_rising[0].airmass - all_stars_rising[-1].airmass
    range_setting = all_stars_setting[0].airmass - all_stars_setting[-1].airmass
    if range_rising > range_setting:
        pick_list = all_stars_rising
        step = range_rising / 5
    else:
        pick_list = all_stars_setting
        step = range_setting / 5

    target = pick_list[0].airmass
    last_delta = step
    last_star = None
    result = []
    for star in pick_list:
        delta = fabs(target - star.airmass)
        if delta > last_delta:
            result.append(last_star)
            target = target - step
            delta = fabs(target - star.airmass)
        last_star = star
        last_delta = delta
    if len(result) < 6:
        result.append(pick_list[-1])

    return result


def write_run(stars: Sequence[ExtStar], out_name: Union[str, bytes, PathLike], filters: Sequence[int] = (2, 3)) -> None:
    out_path = Path(out_name)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    run = peprun.Run()
    for filter_idx in filters:
        run.enable_filter(filter_idx, True)
    for star in stars:
        run.add_item(star.sao, 'EXT')
    with out_path.open('w') as out:
        json.dump(run, out, cls=peprun.RunEncoder)


def report(stars: Sequence[ExtStar]):
    if stars[0].rising:
        print('Stars are RISING.')
    else:
        print('Stars are SETTING.')

    print()
    print('ID        SAO No    B     V     Airmass')
    print('---------------------------------------')
    for star in stars:
        s_id = star.id
        sao = star.sao
        mag_b = star.mag_b
        mag_v = star.mag_v
        airmass = star.airmass
        print(f'{s_id:9s} {sao:9s} {mag_b:5.3f} {mag_v:5.3f} {airmass:5.3f}')


loc = EarthLocation.from_geodetic(Longitude('8d32m12s'), Latitude('50d6m17s'), height=110.0)
dec_lim = 25.0
obs_time = Time.now()
filter_indexes = None
latitude = None
longitude = None
altitude = 0.0
run_name = None

try:
    opts, args = getopt.getopt(sys.argv[1:], 'f:o:a:h:t:r:d:')
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

for opt, arg in opts:
    if opt == '-f':
        if filter_indexes is None:
            filter_indexes = []
        filter_indexes.append(int(arg))
    elif opt == '-o':
        longitude = Longitude(arg)
    elif opt == '-a':
        latitude = Latitude(arg)
    elif opt == '-h':
        altitude = float(arg)
    elif opt == '-t':
        obs_time = Time(arg, scale='utc')
    elif opt == '-r':
        run_name = arg
    elif opt == '-d':
        dec_lim = float(arg)
    else:
        print(f'Unsupported option {opt}.')
        sys.exit(1)

if filter_indexes is None:
    filter_indexes = (2, 3)
if longitude is not None and latitude is not None:
    loc = EarthLocation.from_geodetic(longitude, latitude, altitude)
s, r = get_ra_limits(obs_time, loc, dec_lim)
objs = query_simbad_objects(s, r, dec_lim)
stars_raw = pick_stars(objs, obs_time, loc)
report(stars_raw)

if run_name is not None:
    write_run(stars_raw, run_name, filter_indexes)
