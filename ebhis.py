from .datacube import Datacube
from .mixins import DatacubeMoments, DatacubeSpectra


class EBHISDatacube(Datacube, DatacubeMoments, DatacubeSpectra):

    def __init__(self, *args, **kwargs):

        super(EBHISDatacube, self).__init__(*args, **kwargs)

        self.header['CUNIT3'] = 'm/s'
        self.header['CTYPE3'] = 'VRAD'
        self.header['SPECSYS'] = 'LSRK'
