# -*- python -*-

"""This module contains tools to test forms of equivalence that are particular cases of CCZ-equivalence, including CCZ-equivalence itself.

"""

from .cython_functions import \
    thickness_spectrum, \
    ccz_equivalent_function, \
    enumerate_ea_classes

from .affine_equivalence import *
