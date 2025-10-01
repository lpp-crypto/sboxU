# -*- python -*-


# This file contains the declarations of all the C++ function that we
# want to expose in this module, and of all the cython functions that
# are provided in ./cython_functions.pyx.


# !SECTION! Imports

from sboxUv2.cython_types cimport *
from sboxUv2.core cimport *
from sboxUv2.core.spectrum cimport cpp_Spectrum,Spectrum


# !SECTION! Declaring the C++ Code


# !SUBSECTION! Differential properties

cdef extern from "../cpp/statistics/differential.hpp":
    cpp_Spectrum cpp_differential_spectrum(
        const cpp_S_box & s,
        const int64_t n_threads
    )
    std_vector[std_vector[int64_t]] cpp_ddt(
        const cpp_S_box & s
    )
    bool cpp_is_differential_uniformity_smaller_than(
        const cpp_S_box & s,
        const int64_t u
    )
    std_vector[std_vector[std_vector[BinWord]]] cpp_xddt(
        const cpp_S_box & s
    )


cdef extern from "../cpp/statistics/differential.cpp":
    pass


# !SUBSECTION! Linear properties

cdef extern from "../cpp/statistics/linear.hpp":
    std_vector[int64_t] cpp_walsh_transform(
        const cpp_S_box & s
    )
    cpp_Spectrum cpp_walsh_spectrum(
        const cpp_S_box & s,
        const int64_t n_threads
    )
    std_vector[std_vector[int64_t]] cpp_lat(
        const cpp_S_box & s
    )
    cpp_S_box cpp_invert_lat(
        std_vector[std_vector[int64_t]] & table
    )


cdef extern from "../cpp/statistics/linear.cpp":
    pass



# !SUBSECTION! Boomerang properties

cdef extern from "../cpp/statistics/boomerang.hpp":
    cpp_Spectrum cpp_boomerang_spectrum(
        const cpp_S_box & s,
        const int64_t n_threads
    )
    std_vector[std_vector[int64_t]] cpp_bct(
        const cpp_S_box & s
    )
    cpp_Spectrum cpp_fbct_spectrum(
        const cpp_S_box & s,
        const unsigned int n_threads
    )
    std_vector[std_vector[int64_t]] cpp_fbct(
        const cpp_S_box & s
    )


cdef extern from "../cpp/statistics/boomerang.cpp":
    pass



