# -*- python -*-

from sboxUv2.cython_types cimport *


# !SECTION! Declaring the C++ Code


cdef extern from "../../cpp/core/s_box.hpp":

    # !SUBSECTION! The cpp_S_box class
    
    cppclass cpp_S_box:
        cpp_S_box()
        cpp_S_box(const cpp_S_box s)
        cpp_S_box(std_vector[BinWord] lut)
        cpp_S_box(Bytearray b)
        cpp_S_box(std_vector[BinWord] lut, int64_t input_length, int64_t output_length)
        void destruct()
        BinWord brackets "operator[]" (const BinWord x) const
        cpp_S_box add "operator+" (const cpp_S_box s) const
        cpp_S_box mul "operator*" (const cpp_S_box s) const
        string content_string_repr()
        Bytearray to_bytes() const
        std_vector[BinWord] get_lut()
        int64_t size()
        int64_t get_input_length()
        int64_t input_space_size()
        int64_t get_output_length()
        int64_t output_space_size()
        bool is_invertible()
        cpp_S_box inverse()
        cpp_S_box component(BinWord a)
        cpp_S_box coordinate(BinWord a)
        cpp_S_box derivative(BinWord delta)

    # !SUBSECTION!  S_box generating functions
    
    cpp_S_box cpp_translation(const BinWord a, const int64_t input_bit_length)
    # cpp_S_box cpp_rand_S_box(
    #     cpp_PRNG & alea,
    #     unsigned int input_length,
    #     unsigned int output_length)


# !SUBSECTION! Loading the cpp file 
cdef extern from "../../cpp/core/s_box.cpp":
    pass

# !SUBSECTION! The cpp_S_box_fp class, header only hence no cpp file down below

cdef extern from "../../cpp/core/s_box_fp.hpp":
    cppclass cpp_S_box_fp:
        cpp_S_box_fp()
        cpp_S_box_fp(BinWord input_size, BinWord output_size, cpp_Integer p, std_vector[cpp_Integer] powers_in, std_vector[cpp_Integer] powers_out, std_vector[FpWord] input_space, std_vector[FpWord] output_space, std_vector[FpWord] lut)
        cpp_S_box_fp(cpp_Integer p, std_vector[FpWord] lut)
        BinWord get_input_size() const
        BinWord get_output_size() const
        cpp_Integer get_p() const
        const std_vector[cpp_Integer]& get_powers_in() const
        const std_vector[cpp_Integer]& get_powers_out() const
        const std_vector[FpWord]& get_input_space() const
        const std_vector[FpWord]& get_output_space() const
        const std_vector[FpWord]& get_lut() const

        FpWord operator[](const FpWord& input) const
        cpp_S_box_fp operator+(const cpp_S_box_fp& s) except +  
        cpp_S_box_fp operator*(const cpp_S_box_fp& s) except + 

        bool is_invertible() const

        cpp_S_box_fp get_inverse() const

        cpp_S_box_fp derivative(const FpWord& delta) const

        cpp_S_box_fp coordinate(const BinWord i) const 

        @staticmethod
        std_vector[cpp_Integer] iterated_powers(cpp_Integer p, cpp_Integer n)
        @staticmethod
        std_vector[FpWord] build_input_space(cpp_Integer p, BinWord inputsize)
        @staticmethod 
        FpWord int_to_vec(cpp_Integer i, const std_vector[FpWord]& lookup)
        @staticmethod
        cpp_Integer vec_to_int(const FpWord& v, const std_vector[cpp_Integer]& powers)

# !SECTION! Declaring cython functions and classes

# !SUBSECTION! The S_box class

### TODO : rename as S_box_bin

cdef class S_box:
    cdef cpp_S_box * cpp_sb
    cdef string cpp_name
    cdef list input_cast
    cdef list output_cast
    cdef set_inner_sbox(S_box self, cpp_S_box s)


# !SUBSECTION! Wrapper for the S_box class operators
cdef cpp_S_box pyx_add_sboxes(cpp_S_box s, cpp_S_box t)
cdef cpp_S_box pyx_mul_sboxes(cpp_S_box s, cpp_S_box t)


# !SUBSECTION! The S_Box_fp class

cdef class S_box_fp:
    cdef cpp_S_box_fp * cpp_sb 
    cdef string cpp_name
    cdef set_inner_sbox(S_box_fp self, cpp_S_box_fp s)
