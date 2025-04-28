#!/usr/bin/sage
# Time-stamp: <2025-04-28 16:59:39>

from sage.all import *
from sage.crypto.sbox import SBox
import itertools
from collections import defaultdict

from .sboxU_cython import *
from .linear import *



def get_block_lengths(s):
    """Return the number of bits `n` in the input and `m` in the
    output of `s`

    """
    # finding n
    n = 1
    while (1 << n) < len(s):
        n += 1
    if 2**n != len(s):
        raise Exception("wrong S-box length")
    else:
        # finding m
        mask = 1
        m = 1
        for x in s:
            while (x != (x & mask)):
                mask = (mask << 1) | 1
                m += 1
        return n, m


def preprocess_into_list(s):
    """Returns the content of `s` as a list of `int`.

    The C++ part of sboxU does not understand types such as `SBox` or
    `sage.ring.integer.Integer`. The purpose of this function is to
    get rid of those.

    """
    if isinstance(s, sage.crypto.mq.SBox):
        if type(s[0]) == int:
            return list(s)
        else:
            return [int(x for x in s)]
    else:
        if type(s[0]) == int:
            return s
        else:
            return [int(x for x in s)]
        

def random_function_of_degree(n, m, deg):
    """Returns a function picked randomly in the set of functions mapping
    n bits to m with algebraic degree at most deg.

    """
    result = [0 for x in range(0, 2**n)]
    for output_bit in range(0, m):
        e_i = int(1 << output_bit)
        for monomial_mask in range(1, 2**n):
            if (hamming_weight(monomial_mask) <= deg) and randint(0, 1) == 1:
                for x in range(0, 2**n):
                    if (x & monomial_mask) == monomial_mask:
                        result[x] = oplus(result[x], e_i)
    # r = PolynomialRing(GF(2**n, name="a"), 'x', n)
    # x_is = r.gens()
    # for output_bit in range(0, m):
    #     pol = r.zero()
    #     if randint(0, 1) == 1:
    #         pol = r.one()
    #     for d in range(1, deg+1):
    #         for monomial_x_is in itertools.combinations(x_is, d):
    #             if randint(0, 1) == 1:
    #                 monomial = r.one()
    #                 for x_i in monomial_x_is:
    #                     monomial = x_i*monomial
    #                 pol += monomial
    #     for y in range(0, 2**n):
    #         f_y_bit = pol([(y >> j) & 1 for j in range(0, n)])
    #         result[y] = result[y] | (int(f_y_bit) << output_bit)
    # del(r)
    return result

        
def image(f):
    img = defaultdict(int)
    for x in range(0, len(f)):
        img[f[x]] += 1
    return img.keys()


def all_fields_of_degree(n):
    """Returns a list of all the fields of characteristic 2 and degree
    `n`, meaning that all the primitive polynomials of degree `n` of
    GF(2) are used to create GF instances.

    """
    result = []
    for p in GF(2).polynomial_ring().polynomials(of_degree=n):
        if p.is_primitive():
            result.append(GF(2**n, modulus=p, name="a"))
    return result


def covered_set(mask):
    """Returns a list containing all the integers whose binary
    representation is covered by `mask`, i.e. covered_set(3) = [0, 1, 2,
    3]."""
    if mask == 0:
        return [0]
    elif (mask & 1) == 0:
        return [(x << 1) for x in covered_set(mask >> 1)]
    else:
        c = [(x << 1) for x in covered_set(mask >> 1)]
        return c + [(x | 1) for x in c]


def lg2(x):
    return float(log(x, 2))



# !SECTION! Finite field arithmetic
# =================================


SAGE_VERSION = tuple([int(x) for x in sage.version.version.split(".")])

if SAGE_VERSION < (9, 8):
    def get_convertors(gf):
        if gf.characteristic() == 2:
            return gf.fetch_int, lambda x : x.integer_representation() 
        else:
            return gf.__call__, lambda x : Integer(x)
    
    def ffe_from_int(gf, x):
        if gf.characteristic() > 2:
            return gf(x)
        else:
            return gf.fetch_int(x)

    def ffe_to_int(x):
        if x == 0:         # necessary because of inconsistent casting
            return x
        elif x == 1:       # same...
            return x
        elif x.base_ring().characteristic() > 2:
            return Integer(x)
        else:
            return x.integer_representation()
        
else:
    def get_convertors(gf):
        return gf.fom_integer, lambda x : x.to_integer


    def ffe_from_int(gf, x):
        return gf.from_integer(x)

    def ffe_to_int(x):
        return x.to_integer()



# !SECTION! Convenient function manipulation
# ==========================================

def inverse(s):
    """Returns the functional inverse of the permutation s."""
    result = [0 for i in range(0, len(s))]
    for x in range(0, len(s)):
        result[s[x]] = x
    return result


def eval_function_like(x, o):
    """Allows querying function-like objects such as lists (LUTs),
    matrices, etc. on a given input in a generic way.

    """
    if not isinstance(x, (int, Integer)):
        raise Exception("trying to evaluate a function on a non-integer input\nl={}\nx={}".format(o, x))
    if isinstance(o, (list, sage.crypto.sbox.SBox)):
        return o[x]
    elif isinstance(o, sage.matrix.matrix0.Matrix):
        return apply_bin_mat(x, o)
    elif "base_ring" in dir(o):
        gf = o.base_ring()
        if x > len(gf):
            raise Exception("input ({:d}) is too large to be cast to an element of {}".format(x, gf))
        else:
            return ffe_to_int(o(ffe_from_int(gf, x)))
    elif "__call__" in dir(o):
        return o(x)
    else:
        raise Exception("don't know how to evaluate o={}".format(o))


def get_lut(o, domain_size=None):
    """Returns the lookup table (as a `list`) of the function-like
    object given in the input.

    """
    if isinstance(o, list):
        return o
    elif isinstance(o, sage.matrix.matrix0.Matrix):
        return [apply_bin_mat(x, o) for x in range(0, 2**o.ncols())]
    elif isinstance(o, FastLinearMapping):
        return [o(x) for x in range(0, 2**o.inner_matrix.ncols())]
    elif "base_ring" in dir(o):
        gf = o.base_ring()
        return [ffe_to_int(o(ffe_from_int(gf, x))) for x in range(0, len(gf))]
    elif "__call__" in dir(o):
        if domain_size == None:
            raise Exception("for such ojects ({}), `domain_size` must be specified".format(type(o)))
        else:
            return [o(x) for x in range(0, domain_size)]
    else:
        return Exception("unsupported object type ({})".format(type(o)))
    

def get_input_size(func_list):
    """For a list of "function-like" object, returns the size of the
    set of its inputs (or an exception if a contradiction is found,
    i.e. functions with mismatched input sizes).

    For instance, returns the `len` of a list corresponding to a
    lookup table.

    """
    input_size = None
    # determining function domain
    for f in func_list:
        local_input_size = None
        if isinstance(f, list):
            local_input_size = len(f)
        elif isinstance(f, sage.crypto.sbox.SBox):
            local_input_size = len(list(f))
        elif isinstance(f, FastLinearMapping):
            local_input_size = 2**input_size
        elif isinstance(f, sage.matrix.matrix0.Matrix):
            local_input_size = 2**f.ncols()
        if local_input_size != None:
            if input_size != None and input_size != local_input_size:
                raise Exception("trying to compose functions with mismatched input sizes")
            else:
                input_size = local_input_size
    return input_size

    
def comp(func_list, input_size=None):
    """Implements the composition of functions. Takes as input a list
    of function-like objects, and returns the lookup table of their
    composition.

    The function of highest index is applied first, then the second to
    last, etc. The aim is to mimic the behaviour of the "o" operator,
    i.e. compose([F, G]) returns the lookup table of `F \\circle G`.

    """
    # determining function domain
    if input_size == None:
        input_size = get_input_size(func_list)
        if input_size == None:
            raise Exception("could not determine input size, `input_size` must be specified")
    # composing functions
    result = list(range(0, input_size))
    for f in reversed(func_list):
        result = [eval_function_like(x, f) for x in result]
    return result


def xor_functions(func_list, input_size=None):
    """Implements an F_2 sum (XOR) of functions. Takes as input a list
    of function-like objects, and returns the lookup table of their
    sum.

    """
    # determining function domain
    if input_size == None:
        input_size = get_input_size(func_list)
        if input_size == None:
            raise Exception("could not determine input size, `input_size` must be specified")
    # composing functions
    result = [0 for x in range(0, input_size)]
    for f in func_list:
        result = [oplus(result[x], eval_function_like(x, f))
                  for x in range(0, input_size)]
    return result
    

def F2_trans(cstte):
    """Returns the function taking as input `x` and returning `x
    \\oplus cstte`, i.e. the translation by `cstte`.

    """
    return lambda x : oplus(x, cstte)


def F_mult(gf, cstte):
    """Returns the function taking as in input `x` and returning `x
    \\otimes cstte`, where the multiplication is done in the finite
    field `gf`.

    If `cstte` is a fraction of the form 1/c, then the multiplication
    is by the inverse *in the finite field* of the element
    corresponding to `c`.

    """
    if isinstance(cstte, (int, Integer)):
        elmt = gf.fetch_int(cstte)
        return lambda x : (gf.fetch_int(x) * elmt).integer_representation()
    elif isinstance(cstte, sage.rings.rational.Rational):
        if cstte.numerator() != 1:
            raise Exception("input must be either an integer or a fraction of the form 1/c")
        else:
            elmt = gf.fetch_int(cstte.denominator())**-1
            return lambda x : (gf.fetch_int(x) * elmt).integer_representation()
    elif isinstance(cstte, float):
        rational_form_numerator = int(round((1.0 / cstte)))
        return F_mult(gf, sage.rings.rational.Rational((1, rational_form_numerator)))
    else:
        raise Exception("invalid input type: {}".format(type(cstte)))
        
    



# !SECTION! S-box encoding as bytes
# =================================

# The relevance of the following functions comes from interacting with a database



def pack_to_bytes(s, m):
    """Packs the content of s into a sequence of 8-bit blocks, i.e. a
    `bytearray`. Assumes that all elements of `s` are strictly smaller
    than 2**m.

    """
    if m <= 4:
        result = [0]*floor(len(s) / 2)
        for i in range(0, len(s), 2):
            result[i >> 1] = (s[i] << 4) | s[i+1]
        if (len(s) % 2) == 1: # handling the case of an odd length
            result.append(s[-1])
        return bytearray(result)
    elif m <= 8:
        return bytearray(s)
    else:
        byte_length = ceil(m / 8)
        result = [0] * len(s) * byte_length
        for i in range(0, len(s)):
            x = s[i]
            for j in range(0, byte_length):
                result[i * byte_length + j] = x & 0xFF
                x = x >> 8
        return bytearray(result)

    
def encode_lut(s, m):
    """Returns an array of bytes encoding the full lookup table `s`.

    Since tinySQL doesn't support arrays, we instead store them as
    BLOBs.

    """
    return  bytes([m]) + pack_to_bytes(s, m)
    

def decode_lut(l):
    """Returns a list of integers corresponding to the lut encoded by
    the bytearray l.

    """
    m = l[0]
    b = [int(x) for x in l[1:]]
    if m <= 4:
        result = [0] * len(b)
        for x in range(0, len(result), 2):
            y = b[(x >> 1)]
            result[x]   = y >> 4
            result[x+1] = y & 0xF
    else:
        block = ceil(m / 8)
        result = [0 for x in range(0, int(len(b) / block))]
        for x in range(0, len(result)):
            result[x] = sum(int(b[x*block + j]) << (8*j) for j in range(0, block))
    return result
    
