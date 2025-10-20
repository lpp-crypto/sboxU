# -*- python -*-

from sboxUv2.cython_types cimport *
from sboxUv2.core cimport *
from sboxUv2.statistics cimport *
from sboxUv2.ccz cimport *


# !SECTION! Declaring C++ code

# !SUBSECTION! Invariants

cdef extern from "../cpp/apn/invariants.hpp":
    cpp_S_box cpp_ortho_derivative(
        const cpp_S_box &q
    )

    cpp_Spectrum cpp_sigma_multiplicities(
        const cpp_S_box &f,
        const int64_t k,
        const int64_t n_threads
    )

    
    string cpp_apn_ea_mugshot(
        const cpp_S_box &s,
        const unsigned int n_threads
    )

    
cdef extern from "../cpp/apn/invariants.cpp":
    pass


# !SUBSECTION! Exploring the CCZ class

cdef extern from "../cpp/apn/ccz_class.hpp":
    std_vector[cpp_S_box] cpp_enumerate_ea_classes_quadratic_apn(
        const cpp_S_box &s,
        const unsigned int n_threads
    )

    cpp_S_box cpp_ccz_equivalent_quadratic_function(
        const cpp_S_box & s,
        const unsigned int n_threads
    )



cdef extern from "../cpp/apn/ccz_class.cpp":
    pass
