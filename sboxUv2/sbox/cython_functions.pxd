# -*- python -*-


# This file contains the declarations of all the C++ function that we
# want to expose in this module, and of all the cython functions that
# are provided in ./cython_functions.pyx.


from libcpp.vector cimport vector as cpp_vector # name change imposed by SAGE's own `vector` type
from libcpp.string cimport string
from libcpp cimport bool
from libc.stdint cimport uint64_t, int64_t
from libc.stdint cimport uint64_t, int64_t


cdef extern from "../cpp/s_box.hpp":
    cdef cppclass cpp_S_box:
        cpp_S_box()
        cpp_S_box(cpp_vector[uint64_t] lut)
        cpp_S_box(cpp_vector[uint64_t] lut, int64_t input_length, int64_t output_length)
        uint64_t brackets "operator[]" (const uint64_t x) const
        cpp_S_box add "operator+" (const cpp_S_box s) const
        cpp_S_box mul "operator*" (const cpp_S_box s) const
        string content_string_repr()
        cpp_vector[uint64_t] get_lut()
        int64_t size()
        int64_t get_input_length()
        int64_t input_space_size()
        int64_t get_output_length()
        int64_t output_space_size()
        bool is_invertible()
        cpp_S_box inverse() except +
        
    cpp_S_box cpp_translation(const uint64_t a, const int64_t input_bit_length)

        
cdef extern from "../cpp/s_box.cpp":
    pass


cdef class S_box:
    """doc string doc string"""
    cdef cpp_S_box * cpp_sb
    cdef string cpp_name
    cdef set_inner_sbox(S_box self, cpp_S_box s)

    

cdef cpp_S_box pyx_add_sboxes(cpp_S_box s, cpp_S_box t)
cdef cpp_S_box pyx_mul_sboxes(cpp_S_box s, cpp_S_box t)


cdef S_box pyx_F2_trans(uint64_t k, alphabet)

