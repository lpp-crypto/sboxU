# -*-python-*- 
# Time-stamp: <2024-04-26 10:26:12 leo>
# Modified on 2025-07-30 by baudrin-j

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.string cimport string
from libc.stdint cimport int64_t, uint64_t
from libcpp.utility cimport pair

cdef extern from "sboxu_cpp_utils.cpp":
    pass
    
cdef extern from "sboxu_cpp_utils.hpp" :
    cdef uint64_t oplus_cpp(uint64_t x, uint64_t y)
    cdef unsigned int hamming_weight_cpp(uint64_t x)
    cdef uint64_t parity_cpp(uint64_t x)
    cdef uint64_t scal_prod_cpp(uint64_t x, uint64_t y)
    cdef vector[uint64_t] component_cpp(uint64_t a, vector[uint64_t] f)
    cdef vector[uint64_t] random_permutation_cpp(unsigned int n)
    cdef bint is_permutation_cpp(vector[uint64_t] S)
    cdef vector[uint64_t] inverse_cpp(vector[uint64_t] s)
    cdef int64_t rank_of_vector_set_cpp(vector[uint64_t] l)
    cdef void check_length_cpp(vector[uint64_t] s)


cdef extern from "sboxu_cpp_diff_lin.cpp":
    pass
    
cdef extern from "sboxu_cpp_diff_lin.hpp" :
    cdef map[int64_t, int64_t] differential_spectrum_fast (const vector[uint64_t] s, const unsigned int n)
    cdef vector[ vector[int64_t] ] ddt_cpp(const vector[uint64_t] s)
    cdef vector[int64_t] ddt_row_cpp(const vector[uint64_t] f, const unsigned int a)
    cdef bool is_ddt_row_max_smaller_than_cpp(const vector[uint64_t] s, const int64_t a, const unsigned int u)
    cdef bool is_differential_uniformity_smaller_than_cpp(const vector[uint64_t] s, const int64_t u)
    cdef vector[map[int64_t, int64_t]] c_differential_spectra_cpp(const vector[uint64_t] s, const vector[uint64_t] l_table, const vector[uint64_t] e_table)
    cdef vector[uint64_t] invert_lat_cpp(const vector[vector[int64_t]] l, const unsigned int n)
    cdef vector[int64_t] walsh_spectrum_coord(const vector[uint64_t] f)
    cdef vector[ vector[int64_t] ] lat_cpp(const vector[uint64_t] s)
    cdef vector[uint64_t] lat_zeroes_cpp(const vector[uint64_t] s, const unsigned int n, const unsigned int n_threads)
    cdef vector[uint64_t] projected_lat_zeroes_cpp(const vector[uint64_t] s, const unsigned int n_threads)
    cdef map[int64_t, int64_t] walsh_spectrum_fast_cpp(const vector[uint64_t] s, const unsigned int n_threads)
    cdef vector[uint64_t] ortho_derivative_fast(const vector[uint64_t]& s)
    cdef vector[vector[int64_t]] bct_cpp(const vector[uint64_t] s);
    cdef map[int64_t, int64_t] bct_spectrum_fast_cpp(const vector[uint64_t] s, const unsigned int n_threads)

cdef extern from "sboxu_cpp_equiv.cpp":
    pass

cdef extern from "sboxu_cpp_equiv.hpp" :
    cdef vector[ vector[uint64_t] ] linear_equivalence_cpp(const vector[uint64_t] f, const vector[uint64_t] g, const bool);
    cdef vector[uint64_t] le_class_representative_cpp(const vector[uint64_t] f);

cdef extern from "sboxu_cpp_equiv_approx.cpp":
    pass

cdef extern from "sboxu_cpp_equiv_approx.hpp" :
    cdef vector[ vector[uint64_t] ] linear_equivalence_approx_cpp(const vector[uint64_t] f, const vector[uint64_t] g, const bool a, const unsigned int m);

cdef extern from "sboxu_cpp_ccz.cpp":
    pass

cdef extern from "sboxu_cpp_ccz.hpp" :
    cdef vector[uint64_t] extract_vector_cpp(const vector[uint64_t]& z, const uint64_t a)
    cdef vector[ vector[uint64_t] ] extract_bases_cpp(vector[uint64_t]& z, const int64_t dimension, const int64_t word_length, int64_t n_threads, const string end_condition)
    cdef vector[ vector[uint64_t] ] extract_affine_bases_cpp(vector[uint64_t]& z, const int64_t dimension, const int64_t word_length, int64_t n_threads, const string end_condition)
    cdef map[uint64_t, int64_t] sigma_multiplicities_cpp(const vector[uint64_t] f, const int k, const int n, const int n_threads)
    
cdef extern from "sboxu_cpp_fp.cpp" :
    pass

cdef extern from "sboxu_cpp_fp.hpp" :
    cdef cppclass CppFptFunction:
        CppFptFunction() except +
        CppFptFunction(int64_t, uint64_t, uint64_t, vector[vector[int64_t]], vector[vector[int64_t]]) except +
        int64_t p
        uint64_t t
        uint64_t u
        vector[int64_t] powers 
        vector[vector[int64_t]] indexes
        vector[vector[int64_t]] values
        vector[int64_t] eval(vector[int64_t])
        vector[int64_t] fpt_ddt_row(int64_t)
        vector[ vector[int64_t] ] fpt_ddt()
        map[int64_t,int64_t] fpt_differential_spectrum_fast(const unsigned int)
        
cdef extern from "sboxu_cpp_fp_lat.cpp":
    pass

cdef extern from "sboxu_cpp_fp_lat.hpp":
    cdef vector[ vector[double] ] fpt_lat(const vector[int] s, const int p, const int m, const int num_threads)
    cdef vector[ double ] fpt_lat_column(const vector[int] s, const int p, const int m, const int b)
    cdef vector[ double ] fpt_lat_row(const vector[int] s, const int p, const int m, const int a)
    cdef double fpt_max_lat(const vector[int] s, const int p, const int m, const int num_threads)
    cdef map[double, int] fpt_walsh_spectrum(const vector[int] s, const int p, const int m, const double epsilon, const int num_threads)        

cdef extern from "sboxu_cpp_partition_preserving_linear_mapping.cpp" :
    pass

cdef extern from "sboxu_cpp_partition_preserving_linear_mapping.hpp":
    cdef pair[vector[uint64_t], vector[uint64_t]] cpp_is_linearly_self_equivalent_from_lat(const vector[vector[int64_t]] lat, const string algo, const unsigned int number_of_threads)
    cdef vector[pair[vector[uint64_t], vector[uint64_t]]] cpp_linear_automorphisms_from_lat(const vector[vector[int64_t]] lat, const string algo, const unsigned int number_of_threads)
