# -*- python -*-

from libc.stdint cimport uint64_t
from libcpp.vector cimport vector

from sboxUv2.sbox.cython_functions cimport *
from sboxUv2.statistics.cython_functions cimport *

cdef extern from "../cpp/apn.hpp":
    cpp_S_box cpp_ortho_derivative(const cpp_S_box q)

    
cdef extern from "../cpp/apn.cpp":
    pass

