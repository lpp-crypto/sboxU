#!/usr/bin/sage
# Time-stamp: <2021-08-11 16:41:46 lperrin>

from sage.all import *
import itertools
from collections import defaultdict

from .sboxU_cython import *



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
        

def inverse(s):
    """Returns the functional inverse of the permutation s."""
    result = [0 for i in range(0, len(s))]
    for x in range(0, len(s)):
        result[s[x]] = x
    return result


    
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

