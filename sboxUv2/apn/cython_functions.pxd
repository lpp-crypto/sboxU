# -*- python -*-

from sboxUv2.cython_types cimport *
from sboxUv2.core cimport *
from sboxUv2.statistics cimport *

cdef extern from "../cpp/apn/invariants.hpp":
    cpp_S_box cpp_ortho_derivative(
        const cpp_S_box q
    )

    cpp_Spectrum cpp_sigma_multiplicities(
        const cpp_S_box f,
        const int64_t k,
        const int64_t n_threads
    )
    
cdef extern from "../cpp/apn/invariants.cpp":
    pass

