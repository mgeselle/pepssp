"""Functionality for reading an AAVSO starparm file"""
from astropy.coordinates import SkyCoord, Angle
from astroquery.simbad import Simbad
import astropy.units as u
from csv import DictReader
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Union


@dataclass
class Star:
    usename: str
    auid: Union[None, str]
    pos: SkyCoord
    spec: Union[None, str]
    vmag: float
    bmv: float


@dataclass
class Ensemble:
    pgm: Star
    cmp: Star
    chk: Star


def _tuple_to_hms(h: int, m: int, s: float) -> str:
    return str(h) + 'h' + str(m) + 'm' + str(s) + 's'


def _tuple_to_dms(d: int, m: int, s: float) -> str:
    raw_result = str(d) + 'd' + str(m) + 'm' + str(s) + 's'
    if d > 0:
        return '+' + raw_result
    else:
        return raw_result


def locate_ensemble(parm_file: Union[str, bytes, PathLike], sao_pgm: str) -> Union[None, Ensemble]:
    """Locates the correct line in the starparm file given the SAO number of the program star"""
    parm_path = Path(parm_file)
    if not parm_path.is_file():
        raise FileNotFoundError("starparm file doesn't exist")
    table = Simbad.query_object(sao_pgm)
    pos = None
    for ra, dec in table.iterrows('RA', 'DEC'):
        pos = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    with parm_path.open('r', encoding='utf-8-sig') as in_fd:
        reader = DictReader(in_fd)
        # Stars are considered identical when they are less than 10" apart.
        sep_thresh = Angle('0 0 10 degrees')
        for line in reader:
            ra_str = _tuple_to_hms(int(line['vrah']), int(line['vram']), float(line['vras']))
            de_str = _tuple_to_dms(int(line['vded']), int(line['vdem']), float(line['vdes']))
            cand_pos = SkyCoord(ra_str, de_str, unit=(u.hourangle, u.deg))
            if pos.separation(cand_pos) < sep_thresh:
                usename = line['usename']
                auid = line['auid']
                pgm_vmag = float(line['vvmag'])
                pgm_bmv = float(line['vbmv'])
                pgm_spec = line['vspec']
                pgm = Star(usename, auid, cand_pos, pgm_spec, pgm_vmag, pgm_bmv)
                cmp_name = line['cname']
                cmp_ra_str = _tuple_to_hms(int(line['crah']), int(line['cram']), float(line['cras']))
                cmp_de_str = _tuple_to_dms(int(line['cded']), int(line['cdem']), float(line['cdes']))
                cmp_pos = SkyCoord(cmp_ra_str, cmp_de_str, unit=(u.hourangle, u.deg))
                cmp_vmag = float(line['cvmag'])
                cmp_bmv = float(line['cbmv'])
                cmp = Star(cmp_name, None, cmp_pos, None, cmp_vmag, cmp_bmv)
                chk_name = line['kname']
                chk_ra_str = _tuple_to_hms(int(line['krah']), int(line['kram']), float(line['kras']))
                chk_de_str = _tuple_to_dms(int(line['kded']), int(line['kdem']), float(line['kdes']))
                chk_pos = SkyCoord(chk_ra_str, chk_de_str, unit=(u.hourangle, u.deg))
                chk_vmag = float(line['kvmag'])
                chk_bmv = float(line['deltabmv'])
                chk = Star(chk_name, None, chk_pos, None, chk_vmag, chk_bmv)

                return Ensemble(pgm, cmp, chk)
    return None


if __name__ == '__main__':
    ens = locate_ensemble('/home/mgeselle/Documents/astro/starparm_2021Aug27.csv', 'SAO 8598')
    print(ens)
