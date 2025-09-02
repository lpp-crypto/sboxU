# -*- python -*-

from sboxUv2.cython_types cimport *
from sboxUv2.core cimport *


# !SECTION! Declaring C++ code

# !SUBSECTION! Subspace searches

cdef extern from "../cpp/algorithms/spaceSearch.hpp":

    std_vector[BinWord] cpp_extract_vector(
        const std_vector[BinWord] & z,
        const BinWord a
    )

    std_vector[std_vector[BinWord]] cpp_extract_bases(
        std_vector[BinWord] & z,
        const int64_t dimension,
        int64_t n_threads,
        const string end_condition
    )

    std_vector[std_vector[BinWord]] cpp_extract_affine_bases(
        std_vector[BinWord] & z,
        const int64_t dimension,
        int64_t n_threads,
        const string end_condition
    )

    
cdef extern from "../cpp/algorithms/spaceSearch.cpp":
    pass



# !SUBSECTION! The cpp_BinLinearBasis class
    
cdef extern from "../cpp/algorithms/binLinearBasis.hpp":
    cppclass cpp_BinLinearBasis:
    
        cpp_BinLinearBasis(
            const std_vector[BinWord] & l
        )
        bool add_to_span(
            BinWord x
        )
        bool is_in_span(
            BinWord x
        ) const
        std_vector[BinWord] get_basis() const
        
        int64_t rank() const

        std_vector[BinWord] span() const


        
cdef extern from "../cpp/algorithms/binLinearBasis.cpp":
    pass


# !SECTION! Declaring cython code

cdef class BinLinearBasis:
    cdef cpp_BinLinearBasis * cpp_lb
