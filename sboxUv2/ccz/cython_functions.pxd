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
    

cdef extern from "../cpp/ccz/explore.cpp":
    pass



cdef extern from "../cpp/ccz/partition_preserving_linear_mapping/pplm.cpp" :
    pass

cdef extern from "../cpp/ccz/partition_preserving_linear_mapping/pplm.hpp":
    cdef pair[std_vector[BinWord], std_vector[BinWord]] cpp_is_linearly_self_equivalent_from_lat(const std_vector[std_vector[int64_t]] lat, const string algo, const unsigned int number_of_threads)
    cdef std_vector[pair[std_vector[BinWord], std_vector[BinWord]]] cpp_linear_automorphisms_from_lat(const std_vector[std_vector[int64_t]] lat, const string algo, const unsigned int number_of_threads)