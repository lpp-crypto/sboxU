#!/usr/bin/sage
# Time-stamp: <2023-10-09 10:35:03 lperrin>

from sage.all import *
from sage.crypto.sbox import SBox
import itertools
from collections import defaultdict

from .sboxU_cython import *
from .linear import *


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
    r = PolynomialRing(GF(2**n, name="a"), 'x', n)
    x_is = r.gens()
    for output_bit in range(0, m):
        pol = r.zero()
        if randint(0, 1) == 1:
            pol = r.one()
        for d in range(1, deg+1):
            for monomial_x_is in itertools.combinations(x_is, d):
                if randint(0, 1) == 1:
                    monomial = r.one()
                    for x_i in monomial_x_is:
                        monomial = x_i*monomial
                    pol += monomial
        for y in range(0, 2**n):
            f_y_bit = pol([(y >> j) & 1 for j in range(0, n)])
            result[y] = result[y] | (int(f_y_bit) << output_bit)
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



# !SECTION! Convenient function manipulation
# ------------------------------------------

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
    # !TODO! Polynomial case 
    # elif isinstance(o, sage.rings.polynomial.polynomial_element.Polynomial):
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
    elif "__call__" in dir(o):
        if domain_size == None:
            raise Exception("for such ojects ({}), `domain_size` must be specified".format(type(o)))
        else:
            return [o(x) for x in range(0, domain_size)]
    else:
        return Exception("unsupported object type ({})".format(type(o)))
    
    
def comp(func_list, input_size=None):
    """Implements the composition of functions. Takes as input a list
    of function-like objects, and returns the lookup table of their
    composition.

    The function of highest index is applied first, then the second to
    last, etc. The aim is to mimic the behaviour of the "o" operator,
    i.e. compose([F, G]) returns the lookup table of `F \\circle G`.

    """
    # determining function domain
    f_0 = func_list[0]
    if isinstance(f_0, list):
        input_size = len(f_0)
    elif isinstance(f_0, sage.crypto.sbox.SBox):
        input_size = len(list(f_0))
    elif isinstance(f_0, FastLinearMapping):
        input_size = 2**f_0.inner_matrix.ncols()
    elif isinstance(f_0, sage.matrix.matrix0.Matrix):
        input_size = 2**f_0.ncols()
    elif input_size == None:
        raise Exception("for functions of type {}, `input_size` must be specified".format(type(f_0)))
    # composing functions
    result = list(range(0, input_size))
    for f in reversed(func_list):
        result = [eval_function_like(x, f) for x in result]
    return result


def xor_lut(func_list, input_size=None):
    """Implements an F_2 sum (XOR) of functions. Takes as input a list
    of function-like objects, and returns the lookup table of their
    sum.

    """
    # determining function domain
    f_0 = func_list[0]
    if isinstance(f_0, list):
        input_size = len(f_0)
    elif isinstance(f_0, sage.crypto.sbox.SBox):
        input_size = len(list(f_0))
    elif isinstance(f_0, FastLinearMapping):
        input_size = 2**f_0.inner_matrix.ncols()
    elif isinstance(f_0, sage.matrix.matrix0.Matrix):
        input_size = 2**f_0.ncols()
    elif input_size == None:
        raise Exception("for functions of type {}, `input_size` must be specified".format(type(f_0)))
    # composing functions
    result = [0 for x in range(0, input_size)]
    for f in func_list:
        result = [oplus(result[x], eval_function_like(x, f))
                  for x in range(0, input_size)]
    return result
    

def F2_trans(cstte):
    """Returns the function taking as in input `x` and returning `x
    \\oplus cstte`, i.e. the translation by `cstte`.

    """
    return lambda x : oplus(x, cstte)
