# -*-python-*- 
# Time-stamp: <2021-08-11 14:47:48 lperrin>

from sboxu_cpp cimport *

def linear_equivalence_fast(l0, l1):
	return linear_equivalence_cpp(l0, l1)


def le_class_representative(l0):
	return le_class_representative_cpp(l0)

