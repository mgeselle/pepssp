import ssp
import threading as thr
import datetime


class MeasurementSpec:
    def __init__(self, callback):
        self.callback = callback
        self.exposure_by_filter = {}

    def add_spec(self, filter_index, exposure):
        self.exposure_by_filter[filter_index] = exposure


class Measurement:
    def __init__(self, filter_index, exposure, timestamp):
        self.filter_index = filter_index
        self.exposure = exposure
        self.timestamp = timestamp
        self.counts = []

    def add_count(self, count):
        self.counts.append(count)


class MeasurementRunner:
    def __init__(self, ssp_path, max_integration_time=5999):
        self.__ssp = ssp.SSP(ssp_path, True)
        self.max_integration_time = max_integration_time
        self.__job_event = thr.Event()
        self.__spec = None
        self.__current_filter = 1
        self.__thread = None

    def start(self):
        if self.__thread:
            return
        self.__ssp.open()
        self.__thread = thr.Thread(target=self.__run)
        self.__thread.start()

    def stop(self):
        if not self.__thread:
            return
        self.__job_event.set()
        self.__thread.join()
        self.__thread = None
        self.__ssp.close()

    def submit(self, spec):
        self.__spec = spec
        self.__job_event.set()

    def __run(self):
        while True:
            while not self.__job_event.is_set():
                self.__job_event.wait(1.0)
            self.__job_event.clear()
            if not self.__spec:
                return
            filter_sequence = list(sorted(self.__spec.exposure_by_filter.keys()))

            if len(filter_sequence) > 1 and self.__current_filter == filter_sequence[-1]:
                filter_sequence.reverse()
            measurements = []
            for filter_index in filter_sequence:
                measurements.append(
                    self.__single_measurement(filter_index, self.__spec.exposure_by_filter[filter_index]))
            callback = self.__spec.callback
            self.__spec = None
            callback(measurements)

    def __single_measurement(self, filter_index, integration_time):
        self.__ssp.set_filter(filter_index)
        self.__current_filter = filter_index
        if integration_time == -1:
            integration_time = self.__get_integration_time()
        self.__ssp.set_integ(integration_time)
        measurement = Measurement(filter_index, integration_time, datetime.datetime.utcnow())
        for i in range(3):
            measurement.add_count(self.__ssp.measure())
        return measurement

    def __get_integration_time(self):
        integration_time = 1000
        counts = 0
        while 1 <= integration_time < self.max_integration_time and counts < 5000:
            self.__ssp.set_integ(integration_time)
            counts = self.__ssp.measure()
            if counts > 9999:
                integration_time = int(integration_time / 2)
            elif counts < 5000:
                integration_time = min(self.max_integration_time, int(integration_time * 5000 / counts))
                break
        return integration_time


if __name__ == '__main__':
    def print_measurements(measurements):
        for m in measurements:
            print('Filter {:d}, integration: {:5.2f}s'.format(m.filter_index, m.exposure / 100.0))
            print('Counts: {}'.format(str(m.counts)))

    runner = MeasurementRunner('/dev/ttyUSB0', max_integration_time=1100)
    spc = MeasurementSpec(print_measurements)
    spc.add_spec(2, -1)
    spc.add_spec(3, -1)
    runner.start()
    runner.submit(spc)
