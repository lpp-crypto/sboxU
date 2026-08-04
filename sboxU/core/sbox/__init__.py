"""This module contains the various utilities needed to store and
generate S-boxes.

The idea here is not yet to study S-boxes, only to generate them, and store them in a way that allows calling C++ functions without 
"""



from sboxU.core.sbox.cython_functions import \
    inverse, is_permutation, \
    get_sbox, S_box, S_box_Fp \
