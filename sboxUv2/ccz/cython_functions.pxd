# -*- python -*-

from sboxUv2.cython_types cimport *

from sboxUv2.core cimport *
from sboxUv2.statistics cimport *
from sboxUv2.algorithms cimport *


# !SECTION! Loading C++ code


cdef extern from "../cpp/ccz/zeroes.hpp":
    cpp_Spectrum cpp_thickness_spectrum(
        const cpp_S_box & s,
        const unsigned int n_threads
    )
    
cdef extern from "../cpp/ccz/zeroes.cpp":
    pass


cdef extern from "../cpp/ccz/graph.hpp":
    pass

cdef extern from "../cpp/ccz/graph.cpp":
    pass


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



cdef extern from "../cpp/ccz/partition_preserving_linear_mapping/pplm.cpp" :
    pass

cdef extern from "../cpp/ccz/partition_preserving_linear_mapping/pplm.hpp":
    cdef std_vector[std_vector[cpp_BinLinearMap]] cpp_equivalences_from_lat(
        const std_vector[std_vector[int64_t]] lat1,
        const std_vector[std_vector[int64_t]] lat2,
        const bool single_non_trivial_answer,
        const unsigned int number_of_threads,
        const string equivalence_type)