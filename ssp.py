from serial import Serial, SerialException


class SSPException(Exception):
    """Base class for exceptions related to communicating with the photometer"""
    pass


class SSPTimeoutException(SSPException):
    """Signals a timeout in communicating with the photometer"""
    pass


class SSP:
    """Abstraction for an Optec Inc. SSP3 or SSP5 photometer"""
    port_name = None
    has_filter = False
    __port = None
    filter_index = 1
    gain = 2
    integ = 100

    def __init__(self, port, has_filter=True):
        self.port_name = port
        self.has_filter = has_filter

    def open(self):
        """Opens the serial port and initialises the photometer.

        If the photometer has a an automated filter slider, the slider is sent
        to the home position.
        Gain is set to 10x.
        Integration time is set to 1s.

        Throws an SSPException on error.
        """
        try:
            self.__port = Serial(self.port_name, 19200, timeout=1)
            self.__port.write('SSMODE'.encode('ascii'))
            self.__expect_ack(2)
            if self.has_filter:
                self.__port.write('SHOME.'.encode('ascii'))
                self.__expect_ack(10)
            self.set_gain(self.gain)
            self.set_integ(self.integ)
        except SerialException as e:
            raise SSPException('Error opening serial port', e.strerror)

    def set_gain(self, gain):
        """Sets the gain

        Parameters
        ----------
        gain : int
            Gain setting to use.

            * 1: 100x
            * 2: 10x
            * 3: 1x
        """

        if gain not in (1, 2, 3):
            raise SSPException('illegal gain value: {:d}'.format(gain))
        self.__port.write('SGAIN{:d}'.format(gain).encode('ascii'))
        self.__expect_ack(1)
        self.gain = gain

    def set_integ(self, integ):
        """Set the integration time in hundredths of a second.

        Parameters
        ----------
        integ : int
            integration time in hundredths of a second. Valid range is [0, 5999].
            Note according to the documentation the maximum value is 9999. However, using that
            value on my SSP3 made it take sub-second integrations.
        """
        if integ < 0 or integ > 5999:
            raise SSPException('Illegal integration value: {:d}'.format(integ))
        self.__port.write('SI{:04d}'.format(integ).encode('ascii'))
        self.__expect_ack(1)
        self.integ = integ

    def set_filter(self, filter_index):
        """Selects the filter.

        Parameters
        ----------
        filter_index : int
            Filter index, valid range is [1,6]
        """
        if filter_index < 1 or filter_index > 6:
            raise SSPException('Illegal filter number: {:d}'.format(filter_index))
        self.__port.write('SFILT{:d}'.format(filter_index).encode('ascii'))
        self.__expect_ack(5)
        self.filter_index = filter_index

    def measure(self) -> int:
        """Gets a single reading at the current gain and integration setting."""
        self.__port.write('SCOUNT'.encode('ascii'))
        raw_result = self.__read_ssp(int(self.integ / 100) + 1)
        if len(raw_result) != 7 or raw_result[0:2] != 'C=':
            raise SSPException('Unexpected response from photometer: {}'.format(raw_result))
        return int(raw_result[2:])

    def close(self):
        """Closes the connection to the photometer"""
        self.__port.write('SEND..'.encode('ascii'))
        response = self.__read_ssp(1)
        if response != 'END':
            print('Unexpected response on disconnect: {}'.format(response))
        self.__port.close()
        self.__port = None

    def __read_ssp(self, timeout):
        result = ''
        try:
            for i in range(timeout):
                result = result + str(self.__port.read_until(b'\n\r'), 'utf-8')
                if len(result) > 2 and result[-1] == '\r':
                    return result[0:-2]
            raise SSPTimeoutException()
        except SerialException as e:
            raise SSPException('Error reading from photometer', e.strerror)

    def __expect_ack(self, timeout):
        response = self.__read_ssp(timeout)
        if response != '!':
            raise SSPException('unexpected response from photometer: {}'.format(response))
