# -*- python -*-

from sboxUv2.cython_types cimport *
from sboxUv2.core cimport *
from sboxUv2.algorithms cimport *
from sboxUv2.statistics cimport *


# g++ -o main main.cpp -I /usr/local/include/m4ri -L/usr/local/lib/libm4ri.so.1 -lm4ri -O3
cdef extern from "../../cpp/exploration/sn.hpp":
    
    std_vector[std_vector[cpp_S_box]] cpp_non_trivial_sn (cpp_Integer mode,const cpp_S_box & f, cpp_Integer n_eq, cpp_Integer n_samples)

cdef extern from "../../cpp/exploration/sn.cpp":
    pass