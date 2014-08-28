from .datacube import Datacube
from .mixins import DatacubeMoments, DatacubeSpectra


class CODatacube(Datacube, DatacubeMoments, DatacubeSpectra):

    def __init__(self, *args, **kwargs):

        super(CODatacube, self).__init__(*args, **kwargs)

        self.header['CUNIT3'] = 'm/s'
        self.header['CTYPE3'] = 'VRAD'
        self.header['SPECSYS'] = 'LSRK'
