# cython_functions.pxd

from sboxUv2.cython_types cimport *
from sboxUv2.core cimport *
from sboxUv2.core.spectrum cimport cpp_Spectrum

cdef extern from "../../cpp/core/anf.hpp":
    std_vector[BinWord] cpp_anf_component(const cpp_S_box& f)
    # int64_t cpp_degree_component(const cpp_S_box& f)
    # cpp_Spectrum cpp_monomial_degree_spectrum_component(const cpp_S_box& f)
    cpp_Spectrum cpp_degree_spectrum(const cpp_S_box& f)

cdef extern from "../../cpp/core/anf.cpp":
    pass

