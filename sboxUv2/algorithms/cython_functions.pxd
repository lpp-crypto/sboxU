# -*- python -*-

from sboxUv2.core cimport *
from sboxUv2.cython_types cimport *



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
    
        cpp_BinLinearBasis()
        
        cpp_BinLinearBasis(
            const std_vector[BinWord] & l
        )

        cpp_BinLinearBasis add "operator+"(
            const cpp_BinLinearBasis & L
        ) const

        cpp_BinLinearBasis intersection(
            const cpp_BinLinearBasis & L
        ) const
        
        bool add_to_span(
            BinWord x
        )
        
        bool is_in_span(
            BinWord x
        ) const
        
        std_vector[BinWord] get_basis() const
        
        int64_t rank() const

        std_vector[BinWord] span() const


    std_vector[BinWord] cpp_complete_basis(
        const cpp_BinLinearBasis & basis,
        const unsigned int n
    )

    bool cpp_is_sum_full_rank(
        const cpp_BinLinearBasis & b1,
        const cpp_BinLinearBasis & b2
    )

        
cdef extern from "../cpp/algorithms/binLinearBasis.cpp":
    pass

cdef extern from "../cpp/algorithms/BinLinearBigBasis.hpp":
    cppclass cpp_BinLinearBigBasis:

        cpp_BinLinearBigBasis(unsigned int n)

        cpp_BinLinearBigBasis(
            const std_vector[std_vector[BinWord]] & l, unsigned int n
        )
        
        bool add_to_span(
            std_vector[BinWord] x
        )
        
        bool is_in_span(
            std_vector[BinWord] x
        ) const
        
        std_vector[Bytearray] get_basis() const
        
        int64_t rank() const

        unsigned int size() const

cdef extern from "../cpp/algorithms/BinLinearBigBasis.cpp":
    pass


# !SUBSECTION! The cpp_F2LinearSystem class
    
cdef extern from "../cpp/algorithms/linearSystem.hpp":
    cppclass cpp_F2LinearSystem:
        cpp_F2LinearSystem(const BinWord _n_var, const bool echelonize)
        unsigned int size() const
        BinWord rank() const
        bool add_equation(const std_vector[BinWord] & var_indices)
        void remove_solution(const std_vector[BinWord] & sol)
        std_vector[Bytearray] kernel_as_bytes()
        std_vector[Bytearray] kernel_as_bits()
        string to_string() const

        
cdef extern from "../cpp/algorithms/linearSystem.cpp":
    pass

cdef extern from "../cpp/algorithms/bigvectors.hpp":
    pass


# !SECTION! Declaring cython code

cdef class BinLinearBasis:
    cdef unique_ptr[cpp_BinLinearBasis] cpp_lb


cdef class F2LinearSystem:
    cdef bool echelonize
    cdef unique_ptr[cpp_F2LinearSystem] cpp_ls

cdef class BinLinearBigBasis:
    cdef unique_ptr[cpp_BinLinearBigBasis] cpp_blb
    cdef unsigned int dimension

