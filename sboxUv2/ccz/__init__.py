# -*- python -*-

"""This module contains tools to test forms of equivalence that are particular cases of CCZ-equivalence, including CCZ-equivalence itself.

"""

from .cython_functions import \
    thickness_spectrum, WalshZeroesSpaces, get_WalshZeroesSpaces, \
    ccz_equivalent_function, enumerate_ea_classes, enumerate_permutations_in_ccz_class, EA_mapping, \
    linear_equivalences, el_equivalences, cczl_equivalences, \
    affine_equivalences, ea_equivalences, ccz_equivalences, \
    are_linear_equivalent, are_el_equivalent, are_cczl_equivalent, \
    are_affine_equivalent, are_ea_equivalent, are_ccz_equivalent, \
    ccz_block_decomposition, equivalences_from_lat
from .affine_equivalence import *
from .tu_decompositions import TUdecomposition,thickness,tu_decomposition_from_space_basis,get_tu_decompositions,swap_F2AffineMap
from .equivalences_from_code import are_ccz_equivalent_from_code, are_ea_equivalent_from_code
from .equivalence_from_vq import are_ea_equivalent_from_vq
