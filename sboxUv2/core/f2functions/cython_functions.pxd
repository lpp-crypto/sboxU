# -*- python -*-

from sboxUv2.cython_types cimport *
from ..sbox cimport *


# !SECTION! Declaring C++ code

# !SUBSECTION! Basic functions

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
        const BinWord mask
    )
    int64_t cpp_rank_of_vector_set(
        std_vector[BinWord] l
    )

    std_vector[int] cpp_to_bin ( 
        const BinWord x, 
        int n
    )

    BinWord cpp_from_bin (
        const std_vector[int] & v
    )
    BinWord cpp_circ_shift(
        const BinWord x,
        int n, 
        int shift
    )

    
cdef extern from "../../cpp/core/f2functions.cpp":
    pass


# !SUBSECTION! The cpp_BinLinearMap class

cdef extern from "../../cpp/core/binLinearMap.hpp":
    cppclass cpp_BinLinearMap:
        cpp_BinLinearMap()
    
        cpp_BinLinearMap(
            const std_vector[BinWord] & _image_vectors
        )
    
        cpp_BinLinearMap(
            const std_vector[BinWord] & _image_vectors,
            const int64_t _input_length,
            const int64_t _output_length
        )

        cpp_BinLinearMap(
            const cpp_S_box & lut
        )

        void destruct()

        int64_t get_input_length()
    
        int64_t get_output_length()
    
        BinWord call "operator()" (
            const BinWord x
        ) 
    
        cpp_BinLinearMap mul "operator*"(
           const cpp_BinLinearMap l
        ) 
        
        cpp_BinLinearMap add "operator+"(
            const cpp_BinLinearMap l
        ) 
    
        cpp_BinLinearMap inverse() 
        
        cpp_BinLinearMap transpose()
    
        BinWord rank() 
        
        cpp_S_box get_cpp_S_box() 

        std_vector[BinWord] get_image_vectors()
    
    cpp_BinLinearMap cpp_block_diagonal_BinLinearMap(
        const cpp_BinLinearMap &A,
        const cpp_BinLinearMap &B,
    )



    
cdef extern from "../../cpp/core/binLinearMap.cpp":
    pass


# !SECTION! Declaring cython code


cdef class BinLinearMap:
    cdef cpp_BinLinearMap * cpp_blm

    
