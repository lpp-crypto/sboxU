# -*- python -*-

from sboxU.cython_types cimport *

from sboxU.core cimport *
from sboxU.apn cimport *
from sboxU.statistics cimport *
from sboxU.algorithms cimport *

from libcpp.memory cimport unique_ptr


# !SECTION! Loading C++ code


# !SUBSECTION! Core functionalities

# !SUBSUBSECTION!  zeroes.[hc]pp

cdef extern from "../cpp/ccz/zeroes.hpp":
    cpp_Spectrum cpp_thickness_spectrum(
        const cpp_S_box & s,
        const unsigned int n_threads
    )

    cppclass cpp_WalshZeroesSpaces:
        std_vector[cpp_BinLinearBasis] bases
        
        std_vector[cpp_F2AffineMap] mappings
    
        cpp_WalshZeroesSpaces()

        cpp_WalshZeroesSpaces(
            const std_vector[cpp_BinLinearBasis] & _bases,
            const unsigned int _n,
            const unsigned int _total_size
        )

        cpp_WalshZeroesSpaces(
            const cpp_S_box & s,
            const unsigned int n_threads
        )

        void destruct()

        void init_mappings()

        void init_mappings(const std_vector[cpp_F2AffineMap] & automorphisms)

        void init_mappings(
            const std_vector[cpp_F2AffineMap] & automorphisms_1,
            const std_vector[cpp_F2AffineMap] & automorphisms_2
        )

        cpp_WalshZeroesSpaces image_by(const cpp_F2AffineMap & L) const

        cpp_Spectrum thickness_spectrum() const

    std_vector[std_vector[unsigned int]] cpp_Walsh_zero_orbits(
        const cpp_WalshZeroesSpaces & ws,
        const std_vector[cpp_F2AffineMap] & generators
    )


cdef extern from "../cpp/ccz/zeroes.cpp":
    pass

# !SUBSUBSECTION! graph.[hc]pp

cdef extern from "../cpp/ccz/graph.hpp":
    cdef std_vector[cpp_F2AffineMap] cpp_ccz_block_decomposition(const cpp_F2AffineMap L)

cdef extern from "../cpp/ccz/graph.cpp":
    pass


# !SUBSUBSECTION! explore.[hc]pp

cdef extern from "../cpp/ccz/explore.hpp":
    cpp_S_box cpp_ccz_equivalent_function(
        const cpp_S_box & s,
        const cpp_F2AffineMap & L
        )

    std_vector[cpp_S_box] cpp_enumerate_ea_classes(
        const cpp_S_box & s,
        const unsigned int n_threads
        )

    cpp_F2AffineMap cpp_EA_mapping(
        const cpp_F2AffineMap &A,
        const cpp_F2AffineMap &B,
        const cpp_F2AffineMap &C
        )

    std_vector[cpp_S_box] cpp_enumerate_permutations_in_ccz_class(
        const cpp_S_box & s,
        const unsigned int n_threads
    )

    
cdef extern from "../cpp/ccz/explore.cpp":
    pass


# !SUBSECTION!  Partition preserving mappings

cdef extern from "../cpp/ccz/partition_preserving_linear_mapping/pplm.hpp":
    cdef std_vector[cpp_F2AffineMap] cpp_equivalences_from_lat(
        cpp_S_box sbox1,
        cpp_S_box sbox2,
        const bool single_non_trivial_answer,
        const unsigned int number_of_threads,
        const string equivalence_type)

cdef extern from "../cpp/ccz/partition_preserving_linear_mapping/pplm.cpp" :
    pass

# !SUBSECTION!  Equivalence from Vq

cdef extern from "../cpp/ccz/ea_from_vq.hpp":
    std_vector[cpp_F2AffineMap] cpp_ea_mapping_from_vq(
        const cpp_S_box f,
        const cpp_S_box g,
        const unsigned int n_threads,
        const string & mode)

    std_vector[cpp_F2AffineMap] cpp_product_walsh_match(
        std_vector[cpp_F2AffineMap] G1,
        std_vector[cpp_F2AffineMap] G2,
        cpp_BinLinearBasis V1,
        cpp_BinLinearBasis V2
    )

cdef extern from "../cpp/ccz/ea_from_vq.cpp" :
    pass



# !SECTION! Declaring cython code


# !SUBSECTION! The WalshZeroesSpaces class

cdef class WalshZeroesSpaces:
    cdef unique_ptr[cpp_WalshZeroesSpaces] cpp_wzs;
    cdef list mappings
