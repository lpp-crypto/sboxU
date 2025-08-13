# -*- python -*-

from sboxUv2.cython_types cimport *
from sboxUv2.core cimport *


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

