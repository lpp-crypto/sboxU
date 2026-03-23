# -*- python -*-

"""This module contains tools to investigate APN functions specifically.
"""

from sboxU.apn.cython_functions import \
    ortho_derivative, ortho_integral, \
    sigma_multiplicities, \
    apn_ea_mugshot, apn_ea_mugshot_from_spectra, \
    enumerate_ea_classes_apn_quadratic, ea_mappings_from_ortho_derivative, ccz_equivalent_quadratic_function, \
    automorphisms_from_ortho_derivative, get_WalshZeroesSpaces_quadratic_apn, non_trivial_sn


from sboxU.apn.database import APNFunctions
from sboxU.apn.database import APNQuadraticFunctions_ccz_only


