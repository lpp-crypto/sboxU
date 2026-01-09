# -*- python -*-


from sboxUv2.cython_types cimport *
from sboxUv2.core.sbox cimport *
from sboxUv2.core.spectrum cimport cpp_Spectrum

cdef extern from "../../cpp/core/anf.hpp":
    std_vector[BinWord] cpp_anf_component(const cpp_S_box& f)
    std_vector[BinWord] cpp_quadratic_compact_representation(const cpp_S_box& f)
    std_vector[BinWord] cpp_quadratic_sbox_from_compact_representation(std_vector[BinWord] compact_representation, int64_t n, int64_t m)
    cpp_Spectrum cpp_degree_spectrum(const cpp_S_box &f)
    int64_t cpp_algebraic_degree(const cpp_S_box &f)

cdef extern from "../../cpp/core/anf.cpp":
    pass

