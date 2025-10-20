# -*- python -*-

"""This module contains tools to test forms of equivalence that are particular cases of CCZ-equivalence, including CCZ-equivalence itself.

"""

from .cython_functions import \
    thickness_spectrum, WalshZeroesSpaces, get_WalshZeroesSpaces, \
    ccz_equivalent_function, enumerate_ea_classes, EA_mapping, \
    is_linearly_self_equivalent_from_lat, is_linearly_self_equivalent, \
    linear_automorphisms_from_lat, linear_automorphisms

from .affine_equivalence import *
