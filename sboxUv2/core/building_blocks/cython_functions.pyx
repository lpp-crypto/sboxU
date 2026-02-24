# -*- python -*-

from sage.all import Integer as SAGE_Integer

from cython.operator cimport dereference as ampersand



            
# !SECTION! random generation


# !SUBSECTION! The PRNG itself

cdef class UnsafePRNG:

    def __init__(self, seed : Bytearray=[]):
        if len(seed) == 0:
            seed = cpp_get_seed()
        self.cpp_p = make_unique[cpp_PRNG](<Bytearray>seed)


    def __call__(UnsafePRNG self, vmin: int|SAGE_INTEGER=-1, vmax: int|SAGE_INTEGER=-1) -> BinWord:
        if vmin == -1 and vmax == -1:
            return ampersand(self.cpp_p).call()
        elif vmin == -1 or vmax == -1:
            raise Exception("both vmin and vmax must be set, or none of them")
        else:
            return ampersand(self.cpp_p).call(<BinWord>vmin, <BinWord>vmax)



# !SECTION!  Generating S-Boxes


def rand_invertible_S_box(prng : UnsafePRNG, input_length : int|SAGE_INTEGER) -> S_box:
    result = S_box(name=b"rand_perm")
    (<S_box>result).set_inner_sbox(<cpp_S_box>cpp_rand_invertible_S_box(
        ampersand(prng.cpp_p),
        <int>input_length
    ))
    return result


def rand_S_box(prng : UnsafePRNG, input_length : int|SAGE_INTEGER, output_length : int|SAGE_INTEGER) -> S_box:
    result = S_box(name=b"rand_perm")
    (<S_box>result).set_inner_sbox(<cpp_S_box>cpp_rand_S_box(
        ampersand(prng.cpp_p),
        <int>input_length,
        <int>output_length
    ))
    return result
    
