# -*-python-*- 
# Time-stamp: <2023-01-04 15:57:08 lperrin>

import os
from sboxu_cpp_no_fp_lat cimport *

def linear_equivalence_fast(l0, l1, all_mappings):
    return linear_equivalence_cpp(l0, l1, all_mappings)


def le_class_representative(l0):
    return le_class_representative_cpp(l0)

