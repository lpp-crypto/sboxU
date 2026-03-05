# -*- python -*-

"""Dealing with basic operations over the vector space (F_2)^n (and the finite field F_(2^n).

"""

from sboxUv2.core.f2functions.cython_functions import \
    xor, oplus, \
    hamming_weight, scal_prod, msb, lsb, \
    to_bin, from_bin, circ_shift, \
    linear_combination, rank_of_vector_set, \
    get_F2AffineMap, F2AffineMap, \
    identity_F2AffineMap, zero_F2AffineMap, block_diagonal_F2AffineMap,circ_shift_F2AffineMap

            
from .field_arithmetic import * 

# from .casts import \
#     CastToF2Product, CastFromF2Product
