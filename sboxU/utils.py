#!/usr/bin/sage
# Time-stamp: <2019-05-22 11:37:08 lperrin>

from sage.all import *
from sboxu_cpp import oplus_cpp
import itertools


def oplus(x, y):
    """Ensures that the XOR is computed correctly in a lighter way than
    using Integer(x).__xor__(Integer(y)) as is necessary to work in both
    prompt and script mode.

    """
    return oplus_cpp(long(x), long(y))


def inverse(s):
    """Returns the functional inverse of the permutation s."""
    result = [0 for i in xrange(0, len(s))]
    for x in xrange(0, len(s)):
        result[s[x]] = x
    return result

    
def random_function_of_degree(n, m, deg):
    """Returns a function picked randomly in the set of quadratic
    functions mapping n bits to m.

    """
    result = [0 for x in xrange(0, 2**n)]
    r = PolynomialRing(GF(2**n, name="a"), 'x', n)
    x_is = r.gens()
    for output_bit in xrange(0, m):
        pol = r.zero()
        if randint(0, 1) == 1:
            pol = r.one()
        for d in xrange(1, deg+1):
            for monomial_x_is in itertools.combinations(x_is, d):
                if randint(0, 1) == 1:
                    monomial = r.one()
                    for x_i in monomial_x_is:
                        monomial = x_i*monomial
                    pol += monomial
        for y in xrange(0, 2**n):
            f_y_bit = pol([(y >> j) & 1 for j in xrange(0, n)])
            result[y] = result[y] | (int(f_y_bit) << output_bit)
    return result

        
            
