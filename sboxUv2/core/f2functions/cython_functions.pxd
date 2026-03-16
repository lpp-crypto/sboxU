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


# !SUBSECTION! The cpp_F2AffineMap class

cdef extern from "../../cpp/core/f2affinemap.hpp":
    cppclass cpp_F2AffineMap:
        cpp_F2AffineMap()
    
        cpp_F2AffineMap(
            const std_vector[BinWord] & _image_vectors
        )
    
        cpp_F2AffineMap(
            const std_vector[BinWord] & _image_vectors,
            const int64_t _input_length,
            const int64_t _output_length
        )

        cpp_F2AffineMap(
            const cpp_S_box & lut
        )

        int64_t get_input_length()
    
        int64_t get_output_length()

        bool is_linear()
    
        BinWord operator() (
            const BinWord x
        ) 
    
        cpp_F2AffineMap operator*(
           const cpp_F2AffineMap & l
        ) 
        
        cpp_F2AffineMap operator+(
            const cpp_F2AffineMap & l
        ) 

        bool is_invertible() const
        
        cpp_F2AffineMap inverse() 
        
        cpp_F2AffineMap transpose()
    
        BinWord rank() 
        
        cpp_S_box get_cpp_S_box() 

        std_vector[BinWord] get_image_vectors()
    
    cpp_F2AffineMap cpp_block_diagonal_F2AffineMap(
        const cpp_F2AffineMap &A,
        const cpp_F2AffineMap &B,
    )

    
cdef extern from "../../cpp/core/f2affinemap.cpp":
    pass


# !SECTION! Declaring cython code


cdef class F2AffineMap:
    cdef unique_ptr[cpp_F2AffineMap] cpp_map
    cdef set_inner_map(F2AffineMap self, cpp_F2AffineMap A)

    
