# -*- python -*-

from sboxUv2.cython_types cimport *


# !SECTION! Declaring C++ code

cdef extern from "../../cpp/core/f2functions.hpp":
    uint64_t cpp_msb (
        const uint64_t x
    )
    uint64_t cpp_lsb (
        const uint64_t x
    )
    uint64_t cpp_hamming_weight (
        const uint64_t x
    )
    uint64_t cpp_scal_prod (
        const uint64_t x,
        const uint64_t y
    )
    uint64_t cpp_oplus (
        const uint64_t x,
        const uint64_t y
    )
    uint64_t cpp_linear_combination (
        const cpp_vector[uint64_t] & v,
        uint64_t mask
    )
    int64_t cpp_rank_of_vector_set(
        cpp_vector[uint64_t] l
    )

    # !SUBSECTION! The cpp_Linear_basis class
    
    cppclass cpp_Linear_basis:
    
        cpp_Linear_basis(
            const cpp_vector[uint64_t] & l
        )
        void add_to_span(
            uint64_t x
        )
        bool is_in_span(
            uint64_t x
        ) const
        cpp_vector[uint64_t] get_basis() const
        
        int64_t rank() const

        cpp_vector[uint64_t] span() const

     

    
cdef extern from "../../cpp/core/f2functions.cpp":
    pass


# !SECTION! Declaring cython code


cdef class Linear_basis:
    cdef cpp_Linear_basis * cpp_lb

