# -*- python -*-

from libc.stdint cimport uint64_t
from libcpp.vector cimport vector as cpp_vector

from sboxUv2.sbox.cython_functions cimport *
from sboxUv2.statistics.cython_functions cimport *


cdef extern from "../cpp/algorithms/spaceSearch.hpp":

    cpp_vector[uint64_t] cpp_extract_vector(
        const cpp_vector[uint64_t] & z,
        const uint64_t a
    )

    cpp_vector[cpp_vector[uint64_t]] cpp_extract_bases(
        cpp_vector[uint64_t] & z,
        const int64_t dimension,
        int64_t n_threads,
        const string end_condition
    )

    cpp_vector[cpp_vector[uint64_t]] cpp_extract_affine_bases(
        cpp_vector[uint64_t] & z,
        const int64_t dimension,
        int64_t n_threads,
        const string end_condition
    )

    
cdef extern from "../cpp/algorithms/spaceSearch.cpp":
    pass

