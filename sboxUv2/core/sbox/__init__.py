"""This module contains the various utilities needed to store and
generate S-boxes.

The idea here is not yet to study S-boxes, only to generate them, and store them in a way that allows calling C++ functions without 
"""



from .cython_functions import \
    Sb, S_box, \
    new_sbox_name, \
    F2_trans, identity_S_box

from .misc import \
    random_permutation_S_box, random_function_S_box, \
    F2_mul, monomial, inverse, \
    is_permutation
