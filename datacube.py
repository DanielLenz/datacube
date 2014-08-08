import itertools as it

import numpy as np

from astropy.io import fits
import astropy.wcs as apywcs
from astropy import units as u

from .utilities import brightness_temperature_jybeam


def first_match(keys, d):
    for k in keys:
        if k in d:
            return d[k]
    raise ValueError('No key in mapping.')


class Datacube(object):

    """
    Base class for spectral line data cubes
    """

    _wcs = None
    _resolution = None
    _data = None
    _header = None

    _dtype = None

    _axis_units = None

    _frequencies = None
    _radio_velocities = None
    _optical_velocities = None

    def __init__(self, path=None, data=None, header=None, **kwargs):

        self._dtype = kwargs.get('dtype', np.float32)

        if path is not None:

            for h in fits.open(path):
                if h.is_image:
                    data = h.data
                    header = h.header
                    break

        if (data is not None) and (header is not None):
            
            self._data = u.Quantity(data, u.Unit(header['BUNIT']), self._dtype)
            self._header = fits.Header(header)

        else:
            raise AttributeError(
                "Either path or data and header have to be set.")

        return None

    @property
    def data(self):
        return self._data

    @property
    def header(self):
        return self._header

    @property
    def hdu(self):
        return fits.ImageHDU(self.data.value, self.header)

    @property
    def resolution(self):
        bmaj = first_match(['BMAJ', 'BMAX'], self._header)
        bmin = self._header['BMIN']
        return u.Quantity([bmaj, bmin], u.deg)

    @property
    def wcs(self):
        if self._wcs is None:
            self._wcs = apywcs.WCS(self.header)
        return self._wcs

    @property
    def spec_wcs(self):
        return self.wcs.sub(['spectral'])

    @property
    def cel_wcs(self):
        return self.wcs.sub(['longitude', 'latitude'])

    @property
    def axis_units(self):
        if self._axis_units is None:
            self._axis_units = [u.Unit(s) for s in self.wcs.wcs.cunit]
        return self._axis_units

    @property
    def radio_velocities(self):
        if self._radio_velocities is None:
            rad_eq = u.doppler_radio(self.rest_frequency)
            self._radio_velocities = self.frequencies.to(u.km / u.s, rad_eq)

        return self._radio_velocities

    @property
    def optical_velocities(self):
        if self._optical_velocities is None:
            opt_eq = u.doppler_optical(self.rest_frequency)
            self._optical_velocities = self.frequencies.to(u.km / u.s, opt_eq)

        return self._optical_velocities

    @property
    def frequencies(self):
        if self._frequencies is None:
            specax = self.wcs.wcs.spec

            channels = np.arange(self.data.shape[::-1][specax])

            specc = self.spec_wcs.wcs_pix2world(channels, 0)[0]
            specc *= self.axis_units[specax]

            # Convert radio or optical velocities to frequencies
            if 'VRAD' in self.wcs.wcs.ctype[specax]:
                eq = u.doppler_radio(self.rest_frequency)
            elif 'VOPT' in self.wcs.wcs.ctype[specax]:
                eq = u.doppler_optical(self.rest_frequency)
            elif 'FREQ' in self.wcs.wcs.ctype[specax]:
                eq = []
            else:
                raise AttributeError(
                    'Unsupported spectral type {:s} in header.'.format(self.wcs.wcs.ctype[specax]))

            self._frequencies = specc.to(u.Hz, eq)

        return self._frequencies

    @property
    def rest_frequency(self):
        return self.wcs.wcs.restfrq * u.Hz

    def radio_velocities_to_channels(self, velocities):
        """
        Return the corresponding channels to the given radio velocities.

        Parameters
        ----------
        velocities : array_like
            Velocities in km/s or as astropy.Quantity to search for

        Returns
        -------
        channels : ndarray
            The channels for the given velocities

        true_velocities : ndarray
            The velocities corresponding to the channels
        """

        # If has unit, convert to radio velocity unit
        # else use dimensionless values
        sv = u.Quantity(velocities)
        if sv.unit is not u.dimensionless_unscaled:
            sv = sv.to(self.radio_velocities.unit).value
        else:
            sv = sv.value

        # If spectral axis is ''reversed'', create sorter array
        if self.radio_velocities[0] > self.radio_velocities[-1]:
            sorter = np.arange(self.radio_velocities.size)[::-1]
        else:
            sorter = None

        channels = np.searchsorted(
            self.radio_velocities.value, sv, sorter=sorter)

        true_velocities = self.radio_velocities[channels]

        return channels, true_velocities

    def convert_brightness_unit(self, new_unit):
        """
        Convert the data cube's brightness unit to the
        given new unit. The data are converted in-place
        and the header is updated to reflect the conversion.

        Parameters
        ----------
        new_unit : astropy.Unit
            Compatible unit to convert
            the data to.
        """
        # Hack since self.resolution.prod() is not
        # supported by astropy for whatever reason.
        res = self.resolution
        bunit_eq = brightness_temperature_jybeam(
            np.pi * res[0] * res[1] / 4. / np.log(2), self.rest_frequency)

        self._data = self._data.to(new_unit, bunit_eq)
        self._header['BUNIT'] = self._data.unit.to_string('fits')