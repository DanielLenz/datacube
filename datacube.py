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
    _projection = None
    _data = None
    _header = None

    _dtype = None

    def __init__(self, path=None, data=None, header=None, **kwargs):
        
        self._dtype = kwargs.get('dtype', np.float32)
        
        if path is not None:
            """
            PyFITS open
            """

            self.data, self.header = fits.get_data(path, header=True)

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
        self._header = header

    def _get_header(self):
        return self._header

    header = property(_get_header, _set_header)

    def _get_projection(self):
        if self._projection is None:
            self._projection = apywcs.WCS(self.header)
        return self._projection

    projection = property(_get_projection)



    def _get_velocities(self):
        pass

    velocities = property(_get_velocities)


    def moment(self, vslice=None, cslice=None, kind=[0], mask=None):
        pass




