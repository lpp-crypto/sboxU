# -*- python -*-

from sage.all import Integer as SAGE_INTEGER
from sboxU.core.f2functions import ffe_to_int, circ_shift, i2f_and_f2i
from sboxU.core.f2affinemap import F2AffineMap, get_F2AffineMap
from sboxU.core.sbox import get_sbox

from cython.operator cimport dereference



            
# !SECTION! random generation


# !SUBSECTION! The PRNG itself

cdef class InsecurePRNG:

    def __init__(self, seed : Bytearray=[]):
        if len(seed) == 0:
            seed = cpp_get_seed()
        self.cpp_p = make_unique[cpp_PRNG](<Bytearray>seed)


    def __call__(InsecurePRNG self, vmin: int|SAGE_INTEGER=-1, vmax: int|SAGE_INTEGER=-1) -> BinWord:
        if vmin == -1 and vmax == -1:
            return dereference(self.cpp_p).call()
        elif vmin == -1 or vmax == -1:
            raise Exception("both vmin and vmax must be set, or none of them")
        else:
            return dereference(self.cpp_p).call(<BinWord>vmin, <BinWord>vmax)



# !SECTION!  Generating S-Boxes


# !SUBSECTION! Simple structures


def identity_S_box(length) -> S_box:
    """Returns an S_box instance corresponding to the identity
    function, i.e. the one mapping x to itself.

    """
    return get_sbox(list(range(0, length)))


cdef S_box pyx_F2_trans(BinWord k, n):
    """Wrapper for the `cpp_translation` function. """
    result = S_box(name="Add_{}".format(k))
    result.set_inner_sbox(cpp_translation(k, n))
    return result


def F2_trans(BinWord additive_cstte, field=None, bit_length=None) -> S_box:
    """Returns an S_box containing the lookup table of a simple XOR over a given field extension of F_2.

    If additive_cstte is an integer, then either `field` or `bit_length` must be set. If it is a field element, both `field` and `bit_length` will be ignored.
    
    Args:
        additive_cstte: the constant to add. Can be a field element or an integer. If an integer, then the field used must be specified.
        field: the field in which the multiplication must be made if `additive_cstte` is an integer.
        bit_length: the bit-length to use for both the input and output if `additive_cstte` is an integer.

    Returns:
        An S_box instance
    """
    if isinstance(additive_cstte, (int, SAGE_INTEGER)):
        k = additive_cstte
        if isinstance(bit_length, (int, SAGE_INTEGER)):
            n = bit_length
        elif "degree" in dir(field): # case of a field
            n = field.degree()
        else:
            inputs = {"field": field, "bit_length": bit_length}
            raise Exception("If `additive_cstte` is an integer then either `field` or `bit_length` must be specified, instead, got {}".format(inputs))
    else: # case where the additive constant is a finite field element
        k = ffe_to_int(additive_cstte)
        n = additive_cstte.parent().degree()
    return pyx_F2_trans(k, n)


# !SUBSECTION! Common field operations as S-boxes



def F2_mul(coeff, field=None):
    """Returns an S_box containing the lookup table of a multiplication in an extension of F_2.
    
    Args:
        coeff: the coefficient by which to multiply. Can be a field element or an integer. If an integer, then the field used must be specified.
        field: the field in which the multiplication must be made. If unspecified, the parent field of `coeff` is used.
    
    """
    if isinstance(coeff, (int, SAGE_INTEGER)):
        if field == None:
            raise Exception("If `c` is an integer then the field must be specified!")
        else:
            i2f, f2i = i2f_and_f2i(field)
            c = i2f(coeff)
    else:
        field = coeff.parent()
        i2f, f2i = i2f_and_f2i(field)
        c = coeff
    return get_sbox(
        [
            f2i(i2f(x) * c)
            for x in range(0, field.cardinality())
        ],
        name="Mul_{}".format(f2i(c))
    )

    

def monomial(d, field):
    """Returns an `S_box` containing the LUT of a monomial operating on the given field.

    Args:
        d: the exponent of the monomial (an integer)
        field: a finite field instance assumed to be of characteristic 2.
    """
    assert field.characteristic() == 2
    assert isinstance(d, (int, SAGE_INTEGER))
    i2f, f2i = i2f_and_f2i(field)
    return get_sbox([f2i(i2f(x)**d) for x in range(0, field.cardinality())],
              name="X^{}".format(d))




# !SUBSECTION! Random S-boxes

def rand_invertible_S_box(prng : InsecurePRNG, input_length : int|SAGE_INTEGER) -> S_box:
    result = S_box(name=b"rand_perm")
    (<S_box>result).set_inner_sbox(<cpp_S_box>cpp_rand_invertible_S_box(
        dereference(prng.cpp_p),
        <int>input_length
    ))
    return result


def rand_S_box(prng : InsecurePRNG, input_length : int|SAGE_INTEGER, output_length : int|SAGE_INTEGER) -> S_box:
    result = S_box(name=b"rand_perm")
    (<S_box>result).set_inner_sbox(<cpp_S_box>cpp_rand_S_box(
        dereference(prng.cpp_p),
        <int>input_length,
        <int>output_length
    ))
    return result
    

# !SECTION! Generating F2 affine maps



def identity_F2AffineMap(int64_t n) -> F2AffineMap:
    return get_F2AffineMap([(1 << i) for i in range(0, n)],n,n)


def zero_F2AffineMap(n : BinWord, m : BinWord) -> F2AffineMap:
    return get_F2AffineMap([0 for i in range(0, n)], n, m)


def block_diagonal_F2AffineMap(A, B) -> F2AffineMap:
    Ablm = get_F2AffineMap(A)
    Bblm = get_F2AffineMap(B)
    result = F2AffineMap(None, [], [])
    result.set_inner_map(cpp_block_diagonal_F2AffineMap(
        dereference((<F2AffineMap>Ablm).cpp_map),
        dereference((<F2AffineMap>Bblm).cpp_map),
    ))
    return result


def F2AffineMap_from_blocks(A, B, C, D) -> F2AffineMap:
    Ablm = get_F2AffineMap(A)
    Bblm = get_F2AffineMap(B)
    Cblm = get_F2AffineMap(C)
    Dblm = get_F2AffineMap(D)
    result = F2AffineMap(None, [], [])
    result.set_inner_map(cpp_F2AffineMap_from_blocks(
        dereference((<F2AffineMap>Ablm).cpp_map),
        dereference((<F2AffineMap>Bblm).cpp_map),
        dereference((<F2AffineMap>Cblm).cpp_map),
        dereference((<F2AffineMap>Dblm).cpp_map),
    ))
    return result


def circ_shift_F2AffineMap(int n, int shift) -> F2AffineMap:
    """A circular shift is the operation of rearranging the entries in a vector, either by moving the final entry to the first position, while shifting all other entries to the next position, or by performing the inverse operation. 

    Args : 
        - n : a positive integer
        - shift : a signed integer
    Returns :
        A F2AffineMap object which encodes the circular shift by 'shift' positions. This linear map is an automorphism of (F_2)^n. As for circ_shift, the LSB-first decomposition of a vector x is shifted to the left if shift is positive and to the right otherwise. 
    """
    return get_F2AffineMap([circ_shift(1 <<i,n,shift) for i in range(0, n)], n, n)


def bit_permutation_F2AffineMap(p) -> F2AffineMap:
    """
    A bit permutation is the operation of rearranging the entries in a F2 vector, according to a given permutation. 
    
    Args :
        p : The lut of the bit permutation.
    
    Returns : 
        A BinLinearMap corresponding to bit permutation associated to p.
    """
    return get_F2AffineMap([1 << p[i] for i in range(len(p))])
