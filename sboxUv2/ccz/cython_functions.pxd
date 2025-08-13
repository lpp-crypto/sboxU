# -*- python -*-

from sboxUv2.cython_types cimport *

from sboxUv2.core cimport *
from sboxUv2.statistics cimport *
from sboxUv2.algorithms cimport *


cdef extern from "../cpp/ccz/zeroes.hpp":
    cpp_Spectrum cpp_thickness_spectrum(
        const cpp_S_box & s,
        const unsigned int n_threads
    )
    
cdef extern from "../cpp/ccz/zeroes.cpp":
    pass




cdef extern from "../cpp/ccz/linear_representative.hpp":
    cpp_S_box cpp_le_class_representative(
        const cpp_S_box & f,
        const BinWord fast
    )

cdef extern from "../cpp/ccz/linear_representative.cpp":
    pass
