# -*- python -*-

from libc.stdint cimport uint64_t
from libcpp.vector cimport vector

from sboxUv2.sbox.cython_functions cimport *
from sboxUv2.statistics.cython_functions cimport *

cdef extern from "../cpp/apn/invariants.hpp":
    cpp_S_box cpp_ortho_derivative(
        const cpp_S_box q
    )

    cpp_Spectrum cpp_sigma_multiplicities(
        const cpp_S_box f,
        const int64_t k,
        const int64_t n_threads
    )
    
cdef extern from "../cpp/apn/invariants.cpp":
    pass

