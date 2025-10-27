# -*- python -*-

"""This module contains tools to investigate APN functions specifically.
"""

from .cython_functions import \
    ortho_derivative, sigma_multiplicities, \
    apn_ea_mugshot, apn_ea_mugshot_from_spectra, \
    enumerate_ea_classes_apn_quadratic, ccz_equivalent_quadratic_function, \
    automorphisms_from_ortho_derivative, get_WalshZeroesSpaces_quadratic_apn

from .database import APNFunctions

