from astropy import units as u
from astropy import constants as const

def brightness_temperature_jybeam(pixel_area, beam_area, rest_frequency):
    """
    Defines the conversion between Jy/beam and "brightness temperature"
    """
    beam = beam_area.to(u.sr).value
    pixel = pixel_area.to(u.sr).value
    nu = rest_frequency.to(u.GHz, u.spectral())

    def convert_Jyb_to_K(x_jybm):
        factor = (2 * const.k_B * u.K * nu**2 / const.c**2).to(u.Jy).value
        return (x_jybm / beam / factor)

    def convert_K_to_Jyb(x_K):
        factor = (u.Jy / (2 * const.k_B * nu**2 / const.c**2)).to(u.K).value
        return (x_K * beam / factor)

    def convert_Jyp_to_K(x_jybm):
        factor = (2 * const.k_B * u.K * nu**2 / const.c**2).to(u.Jy).value
        return (x_jybm / pixel / factor)

    def convert_K_to_Jyp(x_K):
        factor = (u.Jy / (2 * const.k_B * nu**2 / const.c**2)).to(u.K).value
        return (x_K * pixel / factor)

    def convert_Jyb_to_Jyp(x_jybm):
        return x_jybm * pixel / beam

    def convert_Jyp_to_Jyb(x_jyp):
        return x_jyp * beam / pixel

    return [
        (u.Jy / u.beam, u.K, convert_Jyb_to_K, convert_K_to_Jyb),
        (u.Jy / u.pixel, u.K, convert_Jyp_to_K, convert_K_to_Jyp),
        (u.Jy / u.pixel, u.Jy / u.beam, convert_Jyp_to_Jyb, convert_Jyb_to_Jyp),
    ]


def first_match(keys, d):
    for k in keys:
        if k in d:
            return d[k]
    raise ValueError('No key in mapping.')