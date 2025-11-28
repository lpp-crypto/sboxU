# -*- python -*-

"""This module contains tools to test forms of equivalence that are particular cases of CCZ-equivalence, including CCZ-equivalence itself.

"""

from .cython_functions import \
    thickness_spectrum, WalshZeroesSpaces, get_WalshZeroesSpaces, \
    ccz_equivalent_function, enumerate_ea_classes, EA_mapping, \
    linear_equivalences_from_lat, extended_linear_equivalences_from_lat, ccz_linear_equivalences_from_lat, \
    linear_equivalences, extended_linear_equivalences, ccz_linear_equivalences, \
    affine_equivalences, extended_affine_equivalences, ccz_equivalences, \
    are_linear_equivalent, are_extended_linear_equivalent, are_ccz_linear_equivalent
from .affine_equivalence import *
from .code_based_equivalence import are_code_equivalent
from .tu_decompositions import TUdecomposition,thickness,tu_decomposition_from_space_basis,get_tu_decompositions,swap_Blm
