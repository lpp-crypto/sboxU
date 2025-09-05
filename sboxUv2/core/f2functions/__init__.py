# -*- python -*-

"""Dealing with basic operations over the vector space (F_2)^n (and the finite field F_(2^n).

"""

from .cython_functions import \
    xor, oplus, \
    hamming_weight, scal_prod, msb, lsb, \
    to_bin, from_bin, \
    linear_combination, rank_of_vector_set, \
    BinLinearMap, Blm, identity_BinLinearMap, zero_BinLinearMap

            
from .field_arithmetic import * 
