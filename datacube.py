from astropy.io import fits
import astropy.wcs as apywcs

class Datacube(object):
    """
    RTFM

    Input
    -----

    path : 

    data : 

    header : 


    Returns
    -------

    """
    _wcs = None
    _data = None
    _header = None

    _dtype = None

    _calculate_velocities = False
    _velocities = None

    def __init__(self, path=None, data=None, header=None, **kwargs):
        
        self._dtype = kwargs.get('dtype', np.float32)
        
        if path is not None:
            """
            PyFITS open
            """

            self.data, self.header = fits.getdata(path, header=True)

        elif (data is not None) and (header is not None):
            """
            Use provided data
            """
            self.data = data
            self.header = header
        
        else:
            raise AttributeError("Either path or data and header have to be set.")
        

    def _set_data(self, data):
        self._data = np.array(data, dtype=self._dtype)

    def _get_data(self):
        return self._data

    data = property(_get_data, _set_data)

    def _set_header(self, header):
        self._calculate_velocities = True
        self._header = header
        self._header['CTYPE3'] = 'VRAD'
        self._header['SPECSYS'] = 'LSRK'


    def _get_header(self):
        return self._header

    header = property(_get_header, _set_header)

    def _get_wcs(self):
        if self._wcs is None:
            self._wcs = apywcs.WCS(self.header)
        return self._wcs

    wcs = property(_get_wcs)

    def _get_velocities(self):
        if (self._velocities is None) or self._calculate_velocities:
            wcs_sub = wcs.sub(['spectral'])

            self._velocities = wcs_sub.wcs_pix2world(np.arange(data.shape[0]), 0) / 1000.
            self._calculate_velocities = False
            
        return self._velocities
        

    velocities = property(_get_velocities)


    def moment(self, vslice=None, cslice=None, kind=[0], mask=None):
        pass




