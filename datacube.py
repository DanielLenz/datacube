import astropy.io as apyio
import astropy.wcs as apywcs

class Datacube(object):

    _projection = None
    _data = None
    _header = None

    _dtype = None

    def __init__(self, path=None, data=None, header=None, **kwargs):

        if path is not None:
            """
            PyFITS open
            """
        else:
            """
            Use provided data
            """

        self._dtype = kwargs.get('dtype', np.float32)


    def _set_data(self, data):
        self._data = np.array(data, dtype=self._dtype)

    def _get_data(self):
        return self._data

    data = property(_get_data, _set_data)



    def _set_header(self, header):
        pass

    def _get_header(self):
        pass

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