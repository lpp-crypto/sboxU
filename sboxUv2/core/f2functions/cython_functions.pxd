# -*- python -*-

from sboxUv2.cython_types cimport *


# !SECTION! Declaring C++ code

cdef extern from "../../cpp/core/f2functions.hpp":
    BinWord cpp_msb (
        const BinWord x
    )
    BinWord cpp_lsb (
        const BinWord x
    )
    BinWord cpp_hamming_weight (
        const BinWord x
    )
    BinWord cpp_scal_prod (
        const BinWord x,
        const BinWord y
    )
    BinWord cpp_oplus (
        const BinWord x,
        const BinWord y
    )
    BinWord cpp_linear_combination (
        const std_vector[BinWord] & v,
        BinWord mask
    )
    int64_t cpp_rank_of_vector_set(
        std_vector[BinWord] l
    )


    
cdef extern from "../../cpp/core/f2functions.cpp":
    pass


# !SECTION! Declaring cython code



