#!/usr/bin/sage
# Time-stamp: <2018-01-25 16:56:54 lperrin>

# from sage.all import RealNumber, RDF, Infinity, exp, log, binomial, factorial, mq
from sage.all import *
from sage.crypto.boolean_function import BooleanFunction

# Loading fast C++ implemented functions
from sboxu_cpp import *


# !SECTION! Wrapping C++ functions for differential/Walsh spectrum 

# Some constants
BIG_SBOX_THRESHOLD = 1024
DEFAULT_N_THREADS  = 16

def walsh_spectrum(s, n_threads=None):
    if n_threads == None:
        if len(s) > BIG_SBOX_THRESHOLD:
            n_threads = DEFAULT_N_THREADS
        else:
            n_threads = 1
    return walsh_spectrum_fast(s, n_threads)


def differential_spectrum(s, n_threads=None):
    if n_threads == None:
        if len(s) > BIG_SBOX_THRESHOLD:
            n_threads = DEFAULT_N_THREADS
        else:
            n_threads = 1
    return differential_spectrum_fast(s, n_threads)
            


# !SECTION! Properties of functions and permutations

# !SUBSECTION! Probability distributions

def lat_coeff_probability_permutation(m, n, c):
    """Returns the probability that a coefficient of the Walsh spectrum of
    a random bijective permutation mapping m bits to n is equal, in
    absolute value, to c.

    If m != n, raises an error.

    """
    if m != n:
        raise "A permutation cannot map {} bits to {}!".format(m, n)
    if c % 4 != 0:
        return 0
    elif c == 0:
        return RealNumber(binomial(2**(n-1), 2**(n-2))**2) / RealNumber(binomial(2**n, 2**(n-1)))
    else:
        c = c/2
        return RealNumber(2) * RealNumber(binomial(2**(n-1), 2**(n-2) + c/2)**2) / RealNumber(binomial(2**n, 2**(n-1)))

    
def lat_coeff_probability_function(m, n, c):
    """Returns the probability that a coefficient of the Walsh spectrum of
    a random function mapping m bits to n is equal, in absolute value, to
    c.

    If m != n, raises an error.

    """
    if m != n:
        raise "m (={}) should be equal to n (={})!".format(m, n)
    if c % 2 != 0:
        return 0
    if c == 0:
        return RDF(2) * RDF(2**(-2**n) * binomial(2**n, 2**(n-1)))
    else:
        c = c/2
        return RDF(4) * RDF(2**(-2**n) * binomial(2**n, 2**(n-1)+c))


def ddt_coeff_probability(m, n, c):
    """Returns the probability that a coefficient of the DDT of a S-Box
    mapping m bits to n is equal to c.

    This probability is identical for random a permutation and a
    random non-bijective function.

    """
    if c % 2 == 1:
        return 0
    else:
        d = c/2
        if m >= 5 and n-m <= m/4:
            k = 2**(m-n-1)
            return RDF(exp(-k) * k**d / factorial(d))
        else:
            return RDF(binomial(2**(m-1), d) * 2**(-n*d) * (1 - 2**-n)**(2**(m-1)-d))

        
def expected_max_ddt(m, n):
    result = RDF(0.0)
    cumul = [0]
    for k in range(1, 15):
        to_add = sum(ddt_coeff_probability(m, n, 2*z) for z in xrange(0, k+1))
        to_add = to_add**RDF((2**m-1) * (2**n-1))
        cumul.append(to_add)
    result = []
    for i in xrange(1, len(cumul)):
        result.append(cumul[i] - cumul[i-1])
    return result


def expected_max_lat(m, n):
    result = RDF(0.0)
    cumul = [0]
    for k in range(1, 50):
        to_add = sum(lat_coeff_probability_permutation(m, n, 2*z) for z in xrange(0, k+1))
        to_add = to_add**RDF((2**m-1) * (2**n-1))
        cumul.append(to_add)
    result = []
    for i in xrange(1, len(cumul)):
        result.append(cumul[i] - cumul[i-1])
    return result

def expected_max_lat_function(m, n):
    result = RDF(0.0)
    cumul = [0]
    for k in range(1, 70):
        to_add = sum(lat_coeff_probability_function(m, n, z) for z in xrange(0, k+1))
        to_add = to_add**RDF((2**m) * (2**n))
        cumul.append(to_add)
    result = []
    for i in xrange(1, len(cumul)):
        result.append(cumul[i] - cumul[i-1])
    return result
    
    
# !SUBSECTION! Aggregated information from the tables

def probability_of_measure(m, n, v_max, occurrences, proba_func):
    """Returns the logarithm in base 2 of the probability that
    $(2^m-1)(2^n-1)$ trials of an experiment yielding output c with
    probability proba_func(m, n, c) will have a result equal to
    `v_max` at most `occurrences` times and be strictly smaller on all
    other trials.

    """
    p_strictly_smaller = sum(RealNumber(proba_func(m, n, i)) for i in xrange(0, v_max))
    p_equal = RealNumber(proba_func(m, n, v_max))
    result = RealNumber(0)
    n_trials = (2**m-1) * (2**n-1)
    for occ in reversed(xrange(0, occurrences+1)):
        added = RealNumber(binomial(n_trials, occ)) * (p_strictly_smaller)**((n_trials - occ)) * p_equal**(occ)
        if abs(added) < Infinity and added != 0 and (result + added > result):
            result += added
    return float(log(result, 2))


def BP_criteria(s):
    """Returns the quantity used by Biryukov and Perrin to generate
    S-Boxes imitating the properties of the S-Box of Skipjack.

    """
    if not isinstance(s, mq.SBox):
        s = mq.SBox(s)
    lat = s.linear_approximation_matrix()
    result = 0
    for i, j in itertools.product(xrange(1, lat.nrows()),
                                  xrange(1, lat.ncols())):
        result += 2**(abs(lat[i, j]))
    return result


# def differential_uniformity(s):
#     if not isinstance(s, mq.SBox):
#         s = mq.SBox(s)
#     if is_permutation(s):
#         return max(differential_coefficients(s).keys())
#     else:
#         diff = differential_coefficients(s)
#         if diff[2**s.m] > 1:
#             return 2**s.m
#         else:
#             coeffs = diff.keys()
#             coeffs.sort()
#             return coeffs[-2]


# def linearity(s):
#     if not isinstance(s, mq.SBox):
#         s = mq.SBox(s)
#     if is_permutation(s):
#         return int(max(linear_coefficients(s).keys()))
#     else:
#         lin = linear_coefficients(s)
#         if lin[2**(s.m-1)] > 1:
#             return int(2**(s.m-1))
#         else:
#             coeffs = lin.keys()
#             coeffs.sort()
#             return int(coeffs[-2])
    


# def differential_branch_number(s):
#     if not isinstance(s, mq.SBox):
#         s = mq.SBox(s)
#     result = 2**s.n
#     for delta_in in xrange(1, 2**s.n):
#         w_in = hamming_weight(delta_in)
#         line = [0 for x in xrange(0, 2**s.m)]
#         for x in xrange(0, 2**s.m):
#             line[oplus(s[x], s[oplus(x,delta_in)])] += 1
#         for delta_out in xrange(0, 2**s.m):
#             if line[delta_out] != 0:
#                 w = w_in + hamming_weight(delta_out)
#                 if w < result:
#                     result = w
#     return result

# def linear_branch_number(s):
#     if not isinstance(s, mq.SBox):
#         s = mq.SBox(s)
#     result = 2**s.n
#     lat = s.linear_approximation_matrix()
#     for a in xrange(1, 2**s.n):
#         w_a = hamming_weight(a)
#         for b in xrange(0, 2**s.m):
#             if lat[a, b] != 0:
#                 w = w_a + hamming_weight(b)
#                 if w < result:
#                     result = w
#     return result

# def probability_of_goodness(s):
#     """Returns the probabilities that the distribution of the coefficients
#     in the DDT and the LAT of a random permutation are at least as
#     good as those of `s`.

#     The criteria used to define "goodness" is the maximum value in the
#     DDT (resp. LAT) and its number of occurrences. Smaller maximum
#     value is better and, given to S-Boxes with the same max value, a
#     smaller number of occurrences is better.

#     """
#     s = mq.SBox(s)
#     result = {}
#     # Looking at the DDT
#     ddt = s.difference_distribution_matrix()
#     dif_measure = max_val_and_occurrences(ddt)
#     result["differential"] = probability_of_measure(
#         len(s),
#         len(s),
#         dif_measure[0],
#         dif_measure[1],
#         ddt_coeff_probability)
#     if True:
#         # Looking at the LAT
#         lat = s.linear_approximation_matrix()
#         lin_measure = max_val_and_occurrences(lat)
#         result["linear"] = probability_of_measure(
#             len(s),
#             len(s),
#             lin_measure[0],
#             lin_measure[1],
#             lat_coeff_probability_permutation)
#         return result
#     else:
#         result["linear"] = "0"
#         return result
    


# !SUBSECTION! Algebraic properies

@parallel
def algebraic_normal_form_coordinate(s, i):
    """Returns the algebraic normal form of the `i`th coordinate of the
    SBox s.

    """
    coordinate = BooleanFunction([(x >> i) & 1 for x in list(s)])
    return coordinate.algebraic_normal_form()


def algebraic_normal_form(s_in):
    """Returns the algebraic normal form of the coordinates of the S-Box
    `s`.

    If we let each of the n output bits of `s` be a Boolean function
    s_i so that $s(x) = s_{n-1}(x) s_{n-2}(x) ... s_{0}(x)$, then
    returns a list l such that l[i] is the ANF of s_i.

    """
    result = {}
    s = mq.SBox(s_in)
    outputs = algebraic_normal_form_coordinate([(s, b) for b in range(0, s.n)])
    for entry in outputs:
        result[entry[0][0][1]] = entry[1]
    return [result[k] for k in range(0, s.n)]


def algebraic_degree(s):
    """Returns the algebraic degree of `s`, i.e. the maximum algebraic
    degree of its coordinates.

    """
    anf = algebraic_normal_form(s)
    result = 0
    for i in xrange(0, len(anf)):
        result = max(result, anf[i].degree())
    return result


# !SECTION! Inverting a LAT

def invert_lat(l):
    """Assuming that l is a valid function LAT, returns this function.

    Only works for functions with the same input and output length."""
    n = int(log(len(l), 2))
    s = [0 for x in xrange(0, 2**n)]
    for i in xrange(0, n):
        b = 1 << i
        for x in xrange(0, 2**n):
            v = int(1/2 - sum(l[a][b]*(-1)**(scal_prod(a, x)) for a in xrange(0, 2**n)) / 2**(n+1))
            if v == 1:
                s[x] = s[x] | (1 << i)
    return s


# !SECTION! Tests

if __name__ == '__main__':
    s = random_permutation(5)
    print s
    l = lat(s)
    s_prime = invert_lat(l)
    print s_prime

