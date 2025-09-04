# cython_functions.pxd

from sboxUv2.cython_types cimport *

cdef extern from "../../cpp/core/s_box.hpp":
    cdef cppclass cpp_S_box:
        cpp_S_box() except +
        int get_input_length() const
        int get_output_length() const
        std_vector[BinWord] get_lut() const
        cpp_S_box component(int u) const

cdef extern from "../../cpp/statistics/spectrum.hpp":
    cdef cppclass cpp_Spectrum:
        cpp_Spectrum() except +
        int64_t size() const
        void incr(int64_t entry)

cdef extern from "../../cpp/core/anf.hpp":
    std_vector[BinWord] cpp_anf_component(const cpp_S_box& f)
    # int64_t cpp_degree_component(const cpp_S_box& f)
    # cpp_Spectrum cpp_monomial_degree_spectrum_component(const cpp_S_box& f)
    cpp_Spectrum cpp_degree_spectrum(const cpp_S_box& f)

cdef extern from "../../cpp/core/anf.cpp":
    pass

