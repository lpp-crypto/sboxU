# -*- python -*-

from sboxUv2.cython_types cimport *

from libcpp.memory cimport unique_ptr


cdef extern from "../../cpp/core/spectrum.hpp":
    cdef cppclass cpp_Spectrum:
        cpp_Spectrum()

        void destruct()
        
        int64_t maximum() const
        
        void incr(
            const int64_t entry
        )

        void incr_by_amount(
            const int64_t entry,
            const int64_t amount
        )

        void incr_by_counting(
            const std_vector[int64_t] & vector_to_count
        )

        int64_t operator[](
            const int64_t key
        ) const

        # renaming is needed because the compiler otherwise mistakes
        # the "keys" method for the one of a dictionary
        std_vector[int64_t] cpp_keys "keys"() const

        int64_t size() const

        string content_string_repr() const

        cpp_Spectrum absolute() const
    
cdef extern from "../../cpp/core/spectrum.cpp":
    pass


# !SUBSECTION! The Spectrum class

cdef class Spectrum:
    cdef string name 
    cdef unique_ptr[cpp_Spectrum] cpp_sp
    cdef set_inner_sp(Spectrum self, cpp_Spectrum sp)
