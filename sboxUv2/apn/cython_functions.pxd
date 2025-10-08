# -*- python -*-

from sboxUv2.cython_types cimport *
from sboxUv2.core cimport *
from sboxUv2.statistics cimport *
from sboxUv2.ccz cimport *


# !SECTION! Declaring C++ code

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

