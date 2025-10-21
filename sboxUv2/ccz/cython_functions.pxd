# -*- python -*-

from sboxUv2.cython_types cimport *

from sboxUv2.core cimport *
from sboxUv2.statistics cimport *
from sboxUv2.algorithms cimport *


# !SECTION! Loading C++ code


# !SUBSECTION! Core functionalities

# !SUBSUBSECTION!  zeroes.[hc]pp

cdef extern from "../cpp/ccz/zeroes.hpp":
    cpp_Spectrum cpp_thickness_spectrum(
        const cpp_S_box & s,
        const unsigned int n_threads
    )

    cppclass cpp_WalshZeroesSpaces:
        cpp_WalshZeroesSpaces()

        cpp_WalshZeroesSpaces(std_vector[cpp_BinLinearBasis] _bases)

        cpp_WalshZeroesSpaces(
            const cpp_S_box & s,
            const unsigned int n_threads
        )

        cpp_WalshZeroesSpaces image_by(const cpp_BinLinearMap & L) const

        cpp_Spectrum thickness_spectrum() const
    
    
cdef extern from "../cpp/ccz/zeroes.cpp":
    pass

# !SUBSUBSECTION! graph.[hc]pp

cdef extern from "../cpp/ccz/graph.hpp":
    pass

cdef extern from "../cpp/ccz/graph.cpp":
    pass


# !SUBSUBSECTION! explore.[hc]pp

cdef extern from "../cpp/ccz/explore.hpp":
    cpp_S_box cpp_ccz_equivalent_function(
        const cpp_S_box & s,
        const cpp_BinLinearMap & L
        )

    std_vector[cpp_S_box] cpp_enumerate_ea_classes(
        const cpp_S_box & s,
        const unsigned int n_threads
        )

    cpp_BinLinearMap cpp_EA_mapping(
        const cpp_BinLinearMap &A,
        const cpp_BinLinearMap &B,
        const cpp_BinLinearMap &C
        )

    
cdef extern from "../cpp/ccz/explore.cpp":
    pass


# !SUBSECTION!  Partition preserving mappings

cdef extern from "../cpp/ccz/partition_preserving_linear_mapping/pplm.hpp":
    cdef pair[std_vector[BinWord], std_vector[BinWord]] cpp_is_linearly_self_equivalent_from_lat(
        const std_vector[std_vector[int64_t]] lat,
        const string algo,
        const unsigned int number_of_threads
    )
    cdef std_vector[pair[std_vector[BinWord], std_vector[BinWord]]] cpp_linear_automorphisms_from_lat(
        const std_vector[std_vector[int64_t]] lat,
        const string algo,
        const unsigned int number_of_threads
    )


cdef extern from "../cpp/ccz/partition_preserving_linear_mapping/pplm.cpp" :
    pass



# !SECTION! Declaring cython code


# !SUBSECTION! The WalshZeroesSpaces class

cdef class WalshZeroesSpaces:
    cdef cpp_WalshZeroesSpaces * cpp_wzs;
