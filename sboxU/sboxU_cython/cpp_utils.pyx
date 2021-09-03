# -*-python-*- 
# Time-stamp: <2021-08-11 16:41:42 lperrin>

from sboxu_cpp cimport *


def oplus(x, y):
    """Ensures that the XOR is computed correctly in a lighter way than
    using Integer(x).__xor__(Integer(y)) as is necessary to work in both
    prompt and script mode.

    """
    return oplus_cpp(x,y)

def hamming_weight(x):
    return hamming_weight_cpp(x)

def parity(x):
    return parity_cpp(x)

def scal_prod(x, y):
    return scal_prod_cpp(x,y)

def component(a, f):
    return component_cpp(a,f)

def random_permutation(n):
    return random_permutation_cpp(n)

def is_permutation(s):
    return is_permutation_cpp(s)
        
def inverse(l):
    return inverse_cpp(l)

def rank_of_vector_set(l):
    return rank_of_vector_set_cpp(l)
