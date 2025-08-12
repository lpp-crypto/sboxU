# -*- python -*-

from libc.stdint cimport uint64_t
from libcpp.vector cimport vector

from sboxUv2.sbox.cython_functions cimport *
from sboxUv2.statistics.cython_functions cimport *
from sboxUv2.algorithms.cython_functions cimport *


cdef extern from "../cpp/ccz/zeroes.hpp":
    cpp_Spectrum cpp_thickness_spectrum(
        const cpp_S_box & s,
        const unsigned int n_threads
    )
    
cdef extern from "../cpp/ccz/zeroes.cpp":
    pass




cdef extern from "../cpp/ccz/linear_representative.hpp":
    cpp_S_box cpp_le_class_representative(
        const cpp_S_box & f
    )

cdef extern from "../cpp/ccz/linear_representative.cpp":
    pass
