# -*-python-*- 
# Time-stamp: <2023-01-04 15:57:08 lperrin>

import os
from sboxU.sboxU_cython.sboxu_cpp_no_fp_lat cimport *

def linear_equivalence_approx_fast(l0, l1, all_mappings, max_contradictions):
    return linear_equivalence_approx_cpp(l0, l1, all_mappings, max_contradictions)


