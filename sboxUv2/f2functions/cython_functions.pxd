# -*- python -*-

from libc.stdint cimport uint64_t
from libcpp.vector cimport vector


cdef extern from "../cpp/f2functions.hpp":
    cdef uint64_t cpp_msb (uint64_t x)
    cdef uint64_t cpp_lsb (uint64_t x)
    cdef uint64_t cpp_hamming_weight (uint64_t x)
    cdef uint64_t cpp_scalar_prod (uint64_t x, uint64_t y)
    cdef uint64_t cpp_oplus (uint64_t x, uint64_t y)
    cdef uint64_t cpp_linear_combination (vector[uint64_t] v, uint64_t mask)

    
cdef extern from "../cpp/f2functions.cpp":
    pass

