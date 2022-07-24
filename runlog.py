"""Functionality for reading peppy run logs"""
import astropy.time as atime
from csv import DictReader
from math import sqrt
from os import PathLike
from pathlib import Path
from statistics import mean, stdev
from typing import Literal, Sequence, Union


class Measurement:
    """Represents a single measurement in all filters, i.e. 3 star and 3 sky deflections per filter"""
    def __init__(self, star: str, star_type: str, time_start: str):
        self.star = star
        self.star_type = star_type
        self.time = atime.Time(time_start, scale='utc')
        self.cnt = dict()
        self.std_cnt = dict()
        self.end_time = None
        self.airmass = None

    def add_counts(self, filt: Literal['U', 'B', 'V', 'R', 'I'], integ_time: int,
                   c1: int, c2: int, c3: int) -> float:
        divisor = integ_time / 1000
        net_cnts = [c1 / divisor, c2 / divisor, c3 / divisor]
        self.cnt[filt] = mean(net_cnts)
        self.std_cnt[filt] = stdev(net_cnts)
        if self.std_cnt[filt] == 0:
            divisor = 0.00001
        else:
            divisor = self.std_cnt[filt]
        return self.cnt[filt] / divisor

    def add_sky_counts(self, sky_time: str, filt: Literal['U', 'B', 'V', 'R', 'I'], integ_time,
                       c1: int, c2: int, c3: int) -> None:
        divisor = integ_time / 1000
        net_cnts = [c1 / divisor, c2 / divisor, c3 / divisor]
        sky_mean = mean(net_cnts)
        sky_stdev = stdev(net_cnts)
        self.cnt[filt] = self.cnt[filt] - sky_mean
        self.std_cnt[filt] = sqrt((self.std_cnt[filt]**2 + sky_stdev**2) / 2)
        self.end_time = atime.Time(sky_time, scale='utc')

    def calc_time(self):
        """Computes the observation time: this is the mean of the time of the
        first star deflection and the last sky deflection."""
        self.time = self.time + ((self.end_time - self.time) / 2)


def read_measurements(input_name: Union[str, bytes, PathLike]) -> Sequence[Measurement]:
    input_path = Path(input_name)
    if not input_path.exists() or not input_path.is_file():
        raise FileNotFoundError("Input file doesn't exist or isn't a file")
    measurements = []
    with open(input_path) as input_file:
        reader = DictReader(input_file)
        last_star = None
        current = None
        line = 1
        for row in reader:
            filt = row['Filter']
            star = row['StarId']
            if star != last_star:
                if current is not None:
                    current.calc_time()
                    measurements.append(current)
                current = Measurement(star, row['StarType'], row['Timestamp'])
                last_star = star
            if row['IsStar'] == 'True':
                snr = current.add_counts(filt, int(row['IntegrationTime']),
                                         int(row['Count1']), int(row['Count2']), int(row['Count3']))
                if snr < 100:
                    print(f'Line {line}: low SNR: {snr:5.1f}')
            else:
                current.add_sky_counts(row['Timestamp'], filt, int(row['IntegrationTime']),
                                       int(row['Count1']), int(row['Count2']), int(row['Count3']))
            line = line + 1
        current.calc_time()
        measurements.append(current)

        return measurements
