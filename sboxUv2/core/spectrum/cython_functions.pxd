# -*- python -*-

from sboxUv2.cython_types cimport *

cdef extern from "../../cpp/core/spectrum.hpp":
    cdef cppclass cpp_Spectrum:
        cpp_Spectrum()
        int64_t maximum() const
        void incr(
            const int64_t entry
        )
        void incr_by_amount(
            const int64_t entry,
            const int64_t amount
        )
        void incr_by_counting(
            const std_vector[int64_t] & vector_to_count
        )
        int64_t brackets "operator[]"(
            const int64_t key
        )
        std_vector[int64_t] keys() const
        int64_t size() const

    
cdef extern from "../../cpp/core/spectrum.cpp":
    pass


# !SUBSECTION! The Spectrum class

cdef class Spectrum:
    cdef string name 
    cdef cpp_Spectrum * cpp_sp
    cdef set_inner_sp(Spectrum self, cpp_Spectrum sp)