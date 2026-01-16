# -*- python -*-

from sage.all import Matrix, GF, Polynomial

from cython.operator cimport dereference as ampersand



            
# !SECTION! random generation

cdef class UnsafePRNG:

    def __init__(self, seed : Bytearray=[]):
        if len(seed) == 0:
            seed = cpp_get_seed()
        self.cpp_p = make_unique[cpp_PRNG](<Bytearray>seed)


    def __call__(UnsafePRNG self, vmin: int=-1, vmax: int=-1) -> BinWord:
        if vmin == -1 and vmax == -1:
            return ampersand(self.cpp_p).call()
        elif vmin == -1 or vmax == -1:
            raise Exception("both vmin and vmax must be set, or none of them")
        else:
            return ampersand(self.cpp_p).call(<BinWord>vmin, <BinWord>vmax)
