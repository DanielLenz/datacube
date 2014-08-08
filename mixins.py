import numpy as np
import itertools as it


class DatacubeMoments(object):

    def moment(self, vslice=None, cslice=None, kind=0, mask=None):
        """
        Make a moment map from either a velocity or channel slice

        Parameters
        ----------
        vslice : flexible
            Velocity slice as either a speed compatible Quantity with two values
            or a two-element array_like with km/s values.

        cslice : array_like
            Channel slice to make the moment of.

        kind : int
            The type of moment to make. Supported are moment 0 and moment 1

        mask : ndarray
            2D (spatial) or 3D mask for moment creation.


        Returns
        -------
        moment : ndarray
            The requested moment of the data cube
        """

        if vslice is not None:
            cslice, _ = self.radio_velocities_to_channels(vslice)
            cslice[-1] += 1

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
                s_velocities = self.radio_velocities[data_slice][:, None, None]

                m = np.nansum(s_data * s_velocities * mask, 0)
                m /= np.nansum(s_data * mask, 0)
                return m
