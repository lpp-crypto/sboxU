# -*- python -*-

from sboxUv2.cython_types cimport *
from ..sbox cimport *


# !SECTION! Declaring C++ code


# !SUBSECTION! The cpp_PRNG class

cdef extern from "../../cpp/core/prng.hpp":
    cppclass cpp_PRNG:
        cpp_PRNG(Bytearray seed)
        BinWord call "operator()" ()
        BinWord call "operator()" (BinWord begin, BinWord end)
        Lut get_permutation(BinWord card)

    Bytearray cpp_get_seed()

cdef extern from "../../cpp/core/prng.cpp":
    pass



# !SUBSECTION! Random Component Generation

cdef extern from "../../cpp/core/componentsGeneration.hpp":
    cdef cpp_S_box cpp_rand_invertible_S_box(cpp_PRNG alea, unsigned int n)
    

cdef extern from "../../cpp/core/componentsGeneration.cpp":
    pass




# !SECTION! Declaring cython code

    
cdef class UnsafePRNG:
    cdef unique_ptr[cpp_PRNG] cpp_p
