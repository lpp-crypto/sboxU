# -*- python -*-


# This file contains the declarations of all the C++ function that we
# want to expose in this module, and of all the cython functions that
# are provided in ./cython_functions.pyx.


# !SECTION! Imports

# !SUBSECTION! Types

from libcpp.vector cimport vector as cpp_vector # name change imposed by SAGE's own `vector` type
from libcpp.string cimport string
from libcpp cimport bool
from libc.stdint cimport uint64_t, int64_t
from libc.stdint cimport uint64_t, int64_t

# !SUBSECTION! Other submodules of sboxU

from sboxUv2.f2functions.cython_functions cimport *
from sboxUv2.sbox.cython_functions cimport *




# !SECTION! Declaring the C++ Code

# !SUBSECTION! The cpp_Spectrum class


cdef extern from "../cpp/statistics/spectrum.hpp":
    cdef cppclass cpp_Spectrum:
        cpp_Spectrum()
        int64_t maximum() const
        void incr(const int64_t entry)
        void incr_by_amount(const int64_t entry, const int64_t amount)
        void incr_by_counting(const cpp_vector[int64_t] vector_to_count)
        int64_t brackets "operator[]"(const int64_t key)
        cpp_vector[int64_t] keys() const
        int64_t size() const

    
cdef extern from "../cpp/statistics/spectrum.cpp":
    pass


# !SUBSECTION! Differential properties

cdef extern from "../cpp/statistics/differential.hpp":
    cpp_Spectrum cpp_differential_spectrum(const cpp_S_box s,
                                           const int64_t n_threads)
    cpp_vector[cpp_vector[int64_t]] cpp_ddt(const cpp_S_box s)


cdef extern from "../cpp/statistics/differential.cpp":
    pass


# !SUBSECTION! Linear properties

cdef extern from "../cpp/statistics/linear.hpp":
    cpp_vector[int64_t] cpp_walsh_transform(const cpp_S_box s)
    cpp_Spectrum cpp_walsh_spectrum(const cpp_S_box s,
                                    const int64_t n_threads)
    cpp_vector[cpp_vector[int64_t]] cpp_lat(const cpp_S_box s)


cdef extern from "../cpp/statistics/linear.cpp":
    pass



# !SUBSECTION! Boomerang properties

cdef extern from "../cpp/statistics/boomerang.hpp":
    cpp_Spectrum cpp_boomerang_spectrum(const cpp_S_box s,
                                        const int64_t n_threads)
    cpp_vector[cpp_vector[int64_t]] cpp_bct(const cpp_S_box s)


cdef extern from "../cpp/statistics/boomerang.cpp":
    pass


# !SECTION! Declaring cython functions and classes

# !SUBSECTION! The Spectrum class

cdef class Spectrum:
    cdef cpp_Spectrum * cpp_sp
    cdef set_inner_sp(Spectrum self, cpp_Spectrum sp)

