# -*- python -*-

from libcpp.vector cimport vector as cpp_vector # name change imposed by SAGE's own `vector` type
from libcpp.string cimport string
from libcpp cimport bool
from libc.stdint cimport uint64_t, int64_t
from libc.stdint cimport uint64_t, int64_t


cdef extern from "../cpp/sbox.hpp":
    cdef cppclass cpp_SBox:
        cpp_SBox()
        cpp_SBox(cpp_vector[uint64_t] lut)
        cpp_SBox(cpp_vector[uint64_t] lut, int64_t input_length, int64_t output_length)
        uint64_t brackets "operator[]" (const uint64_t x) const
        cpp_SBox add "operator+" (const cpp_SBox s) const
        cpp_SBox mul "operator*" (const cpp_SBox s) const
        string content_string_repr()
        cpp_vector[uint64_t] get_lut()
        int64_t size()
        int64_t get_input_length()
        int64_t input_space_size()
        int64_t get_output_length()
        int64_t output_space_size()
        bool is_invertible()
        cpp_SBox inverse() except +
        
    cpp_SBox cpp_translation(const uint64_t a, const int64_t input_bit_length)

        
cdef extern from "../cpp/sbox.cpp":
    pass


cdef class SBox:
    cdef cpp_SBox * cpp_sb
    cdef string cpp_name
    cdef set_inner_sbox(SBox self, cpp_SBox s)

    

cdef cpp_SBox pyx_add_sboxes(cpp_SBox s, cpp_SBox t)
cdef cpp_SBox pyx_mul_sboxes(cpp_SBox s, cpp_SBox t)


cdef SBox pyx_F2_trans(uint64_t k, alphabet)

