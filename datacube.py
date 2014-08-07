import itertools as it

import numpy as np

from astropy.io import fits
import astropy.wcs as apywcs
from astropy import units as u

class Datacube(object):

    """
    Base class for spectral line data cubes
    """

    _wcs = None
    _hdu = None

    _dtype = None

    _axis_units = None

    _spectral_coordinates = None

    def __init__(self, path=None, data=None, header=None, **kwargs):

        self._dtype = kwargs.get('dtype', np.float32)

        if path is not None:

            for h in fits.open(path):
                if h.is_image:
                    self._hdu = h
                    break

        elif (data is not None) and (header is not None):

            self._hdu = fits.ImageHDU(data=data, header=header)

        else:
            raise AttributeError(
                "Either path or data and header have to be set.")

        return None

    
    def _get_data(self):
        if self._hdu.data.dtype != self._dtype:
            self._hdu.data = self._hdu.data.astype(self._dtype)
        return self._hdu.data

    data = property(_get_data)

    
    def _get_header(self):
        return self._hdu.header

    header = property(_get_header)

    
    def _get_hdu(self):
        return self._hdu

    hdu = property(_get_hdu)

    
    def _get_wcs(self):
        if self._wcs is None:
            self._wcs = apywcs.WCS(self.header)
        return self._wcs

    wcs = property(_get_wcs)

    
    def _get_spec_wcs(self):
        return self.wcs.sub(['spectral'])

    spec_wcs = property(_get_spec_wcs)


    def _get_cel_wcs(self):
        return self.wcs.sub(['longitude', 'latitude'])

    cel_wcs = property(_get_cel_wcs)


    def _get_axis_units(self):
        if self._axis_units is None:
            self._axis_units = [u.Unit(s) for s in self.wcs.wcs.cunit]
        return self._axis_units

    axis_units = property(_get_axis_units)


    def _get_radio_velocities(self):
        specc = self._get_spectral_coordinates()
        rad_eq = u.doppler_radio(self.wcs.wcs.restfrq * u.Hz)
        return specc.to(u.km / u.s, equivalencies=rad_eq)

    radio_velocities = property(_get_radio_velocities)


    def _get_spectral_coordinates(self):
        if self._spectral_coordinates is None:
            specax = self.wcs.wcs.spec
            dataax = self.wcs.wcs.naxis - specax - 1
            channels = np.arange(self.data.shape[dataax])
            self._spectral_coordinates = self.spec_wcs.wcs_pix2world(channels, 0)[0]
            self._spectral_coordinates *= self.axis_units[specax]
        
        return self._spectral_coordinates

    
    def moment(self, vslice=None, cslice=None, kind=0, mask=None):

        if vslice is not None:
            cslice = self.spec_wcs.wcs_world2pix(vslice, 0)[-1]

        if cslice is not None:

            cslice = [int(f(c))
                      for f, c in it.izip([np.floor, np.ceil], cslice)]
            data_slice = slice(*cslice)

            if mask is None:
                mask = 1.
            elif mask.shape == self.data.shape[1:]:
                mask = mask[None]
            elif mask.shape == self.data.shape:
                mask = mask[data_slice]

            s_data = self.data[data_slice]

            if kind == 0:
                return np.nansum(s_data * mask, 0)

            if kind == 1:
                s_velocities = self.velocities[data_slice][:, None, None]

                m = np.nansum(s_data * s_velocities * mask, 0)
                m /= np.nansum(s_data * mask, 0)
                return m


class EBHISDatacube(Datacube):


    def __init__(self, *args, **kwargs):

        super(EBHISDatacube, self).__init__(*args, **kwargs)

        self._hdu.header['CUNIT3'] = 'm/s'
        self._hdu.header['CTYPE3'] = 'VRAD'
        self._hdu.header['SPECSYS'] = 'LSRK'
