# -*- python -*-

from libc.stdint cimport uint64_t, int64_t
from libcpp.vector cimport vector as cpp_vector


# !SECTION! Declaring C++ code

cdef extern from "../cpp/f2functions.hpp":
    uint64_t cpp_msb (
        uint64_t x
    )
    uint64_t cpp_lsb (
        uint64_t x
    )
    uint64_t cpp_hamming_weight (
        uint64_t x
    )
    uint64_t cpp_scal_prod (
        uint64_t x,
        uint64_t y
    )
    uint64_t cpp_oplus (
        uint64_t x,
        uint64_t y
    )
    uint64_t cpp_linear_combination (
        cpp_vector[uint64_t] v,
        uint64_t mask
    )
    int64_t cpp_rank_of_vector_set(
        cpp_vector[uint64_t] l
    )

    # !SUBSECTION! The cpp_Linear_basis class
    
    cppclass cpp_Linear_basis:
    
        cpp_Linear_basis(
            cpp_vector[uint64_t] l
        )
        void add_to_span(
            const uint64_t x
        )
        cpp_vector[uint64_t] get_basis() const
        
        int64_t rank() const

        cpp_vector[uint64_t] span() const

     

    
cdef extern from "../cpp/f2functions.cpp":
    pass


# !SECTION! Declaring cython code


cdef class Linear_basis:
    cdef cpp_Linear_basis * cpp_lb

