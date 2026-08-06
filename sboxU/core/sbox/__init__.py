"""This module contains the various utilities needed to store and
generate S-boxes.

The idea here is not yet to study S-boxes, only to generate them, and store them in a way that allows calling C++ functions without 
"""



from sboxU.core.sbox.cython_functions import \
    get_sbox, S_box, \
    new_sbox_name, \
    F2_trans, identity_S_box

from sboxU.core.sbox.misc import \
    random_permutation_S_box, random_function_S_box, \
    F2_mul, monomial, inverse, \
    is_permutation

from sboxU.core.sbox.linearCasts import \
    CastToF2Product, CastFromF2Product, canonical_cast, \
    CastFromF2n, CastToF2n, casts_from_field
