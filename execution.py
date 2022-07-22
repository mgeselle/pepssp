import measurement
import csv


class RunExecution:
    """Handles the execution of a run."""
    def __init__(self, run, filter_slots, device, callback):
        self.__run = run
        self.__filter_slots = filter_slots
        self.__runner = measurement.MeasurementRunner(device)
        self.__callback = callback
        self.__output_fp = None
        self.__output_writer = None
        self.__item_index = -1
        self.__last_measurements = None
        self.__problem_filters = []
        self.__good_measurements = []

    def start(self, output_path):
        self.__output_fp = open(output_path, mode='w', buffering=1)
        self.__output_writer = csv.writer(self.__output_fp)
        header = ['Timestamp', 'Index', 'StarId', 'StarType', 'IsStar',
                  'Filter', 'IntegrationTime', 'Count1', 'Count2', 'Count3']
        self.__output_writer.writerow(header)
        self.__runner.start()
        self.__next_item()

    def finish(self):
        self.__runner.stop()
        self.__output_fp.close()
        self.__output_writer = None

    def operator_ack(self):
        spec = measurement.MeasurementSpec(self.__measurement_done)
        if self.__last_measurements:
            for m in self.__last_measurements:
                spec.add_spec(m.filter_index, m.exposure)
        elif self.__problem_filters:
            for f in self.__problem_filters:
                # ToDo: put integration time selection strategy into GUI
                # Use 20s for blue, 10s for anything else in the meantime
                spec.add_spec(f, 2000 if self.__filter_slots[f - 1] == 'B' else 1000)
        else:
            for f in self.__run.get_filters():
                spec.add_spec(f, 2000 if self.__filter_slots[f - 1] == 'B' else 1000)
        self.__runner.submit(spec)

    def __next_item(self):
        self.__item_index = self.__item_index + 1
        self.__last_measurements = None
        self.__good_measurements = []
        if self.__item_index == self.__run.len():
            self.__callback('Run complete', False, self.__item_index)
            return

        self.__callback('Go to {}'.format(self.__run.get_item(self.__item_index).star_id), True, self.__item_index)

    def __measurement_done(self, measurements):
        # Try to detect issues due to star wandering outside diaphragm...
        self.__problem_filters = []
        if not self.__last_measurements:
            for m in measurements:
                # Arbitrarily using 5% difference between max and min as criterion.
                # We are aiming for SNR > 100. If individual measurements fluctuate by more than 10%,
                # that goal is unattainable. Photometry with accuracy worse than 0.1mag is useless.
                if max(m.counts) / min(m.counts) > 1.05:
                    self.__problem_filters.append(m.filter_index)
                else:
                    self.__good_measurements.append(m)
        else:
            self.__good_measurements = measurements

        if self.__problem_filters:
            self.__callback('max(count)/min(count) > 1.1: check centering', True, self.__item_index)
            return

        for m in self.__good_measurements:
            row = [m.timestamp.isoformat(timespec='seconds'),
                   str(self.__item_index),
                   self.__run.get_item(self.__item_index).star_id,
                   self.__run.get_item(self.__item_index).star_type,
                   str(not bool(self.__last_measurements)),
                   self.__filter_slots[m.filter_index - 1],
                   m.exposure]
            row = row + m.counts
            self.__output_writer.writerow(row)
        if self.__last_measurements:
            self.__next_item()
        else:
            self.__last_measurements = self.__good_measurements
            self.__good_measurements = []
            self.__callback('Go to sky', True, self.__item_index)
