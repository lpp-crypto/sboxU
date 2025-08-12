"""This module contains the various utilities needed to store and
generate S-boxes.

The idea here is not yet to study S-boxes, only to generate them, and store them in a way that allows calling C++ functions without 
"""

import sys



from .cython_functions import Sb, new_sbox_name, S_box, F2_trans, identity_S_box

from .misc import *
