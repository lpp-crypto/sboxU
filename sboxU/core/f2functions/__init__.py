# -*- python -*-

"""Dealing with basic operations over the vector space (F_2)^n (and the finite field F_(2^n).

"""

            
from sboxU.core.f2functions.field_arithmetic import i2f_and_f2i, ffe_from_int, ffe_to_int

from sboxU.core.f2functions.linearcasts import \
    loop_over_structure, canonical_cast, \
    CastFromF2Product, CastToF2Product, \
    CastFromF2n, CastToF2n, casts_from_field

from sboxU.core.f2functions.cython_functions import \
    xor, oplus, \
    hamming_weight, scal_prod, msb, lsb, \
    to_bin, from_bin, circ_shift, \
    linear_combination, rank_of_vector_set, \
    F2Transformation
