from .datacube import Datacube
from .mixins import DatacubeMoments


class EBHISDatacube(Datacube, DatacubeMoments):

    def __init__(self, *args, **kwargs):

        super(EBHISDatacube, self).__init__(*args, **kwargs)

        self._hdu.header['CUNIT3'] = 'm/s'
        self._hdu.header['CTYPE3'] = 'VRAD'
        self._hdu.header['SPECSYS'] = 'LSRK'
