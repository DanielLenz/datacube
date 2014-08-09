from astropy import units as u
from astropy import constants as const

def brightness_temperature_jybeam(beam_area, rest_frequency):
    """
    Defines the conversion between Jy/beam and "brightness temperature"
    """
    beam = beam_area.to(u.sr).value
    nu = rest_frequency.to(u.GHz, u.spectral())

    def convert_Jyb_to_K(x_jybm):
        factor = (2 * const.k_B * u.K * nu**2 / const.c**2).to(u.Jy).value
        return (x_jybm / beam / factor)

    def convert_K_to_Jyb(x_K):
        factor = (u.Jy / (2 * const.k_B * nu**2 / const.c**2)).to(u.K).value
        return (x_K * beam / factor)

    return [(u.Jy / u.beam, u.K, convert_Jyb_to_K, convert_K_to_Jyb)]


def first_match(keys, d):
    for k in keys:
        if k in d:
            return d[k]
    raise ValueError('No key in mapping.')