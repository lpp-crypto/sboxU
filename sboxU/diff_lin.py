#!/usr/bin/sage
# Time-stamp: <2023-10-09 14:02:19 lperrin>


# from sage.all import RealNumber, RDF, Infinity, exp, log, binomial, factorial,
from sage.all import *
from sage.crypto.boolean_function import BooleanFunction
import itertools
from collections import defaultdict

# Loading fast C++ implemented functions
from .sboxU_cython import *
from .utils import *

# !SECTION! Wrapping C++ functions for differential/Walsh spectrum 

# Some constants
BIG_SBOX_THRESHOLD = 128
DEFAULT_N_THREADS  = 2
DEFAULT_HIGH_PRECISION = 40


def lat_zeroes(s, n_threads=None):
    if n_threads == None:
        if len(s) > BIG_SBOX_THRESHOLD:
            n_threads = DEFAULT_N_THREADS
        else:
            n_threads = 1
    n = int(log(len(s), 2))
    return lat_zeroes_fast(s, n, n_threads)

def proj_lat_zeroes(s, n_threads=None):
    if n_threads == None:
        if len(s) > BIG_SBOX_THRESHOLD:
            n_threads = DEFAULT_N_THREADS
        else:
            n_threads = 1
    return projected_lat_zeroes_fast(s, n_threads)

def boomerang_uniformity(s):
    b = bct(s)
    boom_unif = 0
    for i in range(1, len(s)):
        for j in range(1, len(s)):
            boom_unif = max(boom_unif, b[i][j])
    return boom_unif


def dlct(s):
    n = int(log(len(s), 2))
    table = [[0 for b in range(0, 2**n)] for a in range(0, 2**n)]
    for delta in range(0, 2**n):
        for l in range(0, 2**n):
            table[delta][l] = sum((-1)**(scal_prod(l, oplus(s[x], s[oplus(x, delta)])))
                              for x in range(0, 2**n))
    return table


def c_differential_uniformity_spectrum(s, F, l_table=None, e_table=None):
    """Returns a dictionnary `d` such that `d[u]` is a list containing all
    values `c` such that c != 0, and such that s[x ^ a] ^ c s[x] = b
    has at most u solutions x, the maximum being taken over all pairs
    (a, b).

    """
    c_spectra = c_differential_spectra(s, F, l_table=None, e_table=None)
    result = defaultdict(list)
    for c in range(1, len(s)):
        uniformity = max(c_spectra[c].keys())
        result[uniformity].append(c)
    return result


# !SECTION! Linear structures

def linear_structures(f):
    """Assuming that `f` is the LUT of a Boolean function, returns a pair of
    lists `[l_0, l_1]` such that, for all `a` in l_e (with e in [0,1]),

    f(x+a) + f(x) = e,

    for all x (where `+` corresponds to a XOR). Obviously, 0 is always
    in `l_0`.

    """
    result = [[0], []]
    for a in range(1, len(f)):
        offset = oplus(f[a], f[0])
        valid = True
        for x in range(1, len(f)):
            if oplus(f[oplus(x, a)], f[x]) != offset:
                valid = False
                break
        if valid:
            result[offset].append(a)    
    return result


def linear_structures_vectorial(s):
    """Assuming that `s` is the LUT of a vectorial Boolean function,
    returns a dictionnary `d`, where `d[c]` is pair of lists `[l_0,
    l_1]` such that, for all `a` in l_e (with e in [0,1]),

    c. (s(x+a) + s(x)) = e,

    for all x (where `+` corresponds to a XOR), where `.` is the
    scalar product.

    The dictionnary only contains keys where `l_0` or `l_1` are
    non-trivial, i.e. where `l_0` contains more than just 0 and/or
    where `l_1` is non-empty.

    """
    result = {}
    for c in range(1, len(s)):
        f = component(c, s)
        l = linear_structures(f)
        if len(l[0]) > 1 or len(l[1]) > 0:
            result[c] = l
    return result
    

        
# !SECTION! Properties of functions and permutations

# !SUBSECTION! Probability distributions

def lat_coeff_probability_permutation(m, n, c, precision=DEFAULT_HIGH_PRECISION):

    """Returns the probability that a coefficient of the Walsh spectrum of
    a random bijective permutation mapping m bits to n is equal, in
    absolute value, to c.

    If m != n, raises an error.

    """
    big_precision = RealField(precision)
    if m != n:
        raise "A permutation cannot map {} bits to {}!".format(m, n)
    if c % 4 != 0:
        return 0
    elif c == 0:
        return big_precision(Integer(binomial(2**(n-1), 2**(n-2)))**2) / Integer(binomial(2**n, 2**(n-1)))
    else:
        c = c/2
        return 2 * big_precision(Integer(binomial(2**(n-1), 2**(n-2) + c/2))**2) / Integer(binomial(2**n, 2**(n-1)))

    
def lat_coeff_probability_function(m, n, c, precision=DEFAULT_HIGH_PRECISION):
    """Returns the probability that a coefficient of the Walsh spectrum of
    a random function mapping m bits to n is equal, in absolute value, to
    c.

    If m != n, raises an error.

    """
    big_precision = RealField(precision)
    if m != n:
        raise "m (={}) should be equal to n (={})!".format(m, n)
    if c % 4 != 0:
        return 0
    if c == 0:
        return 2 * big_precision(2**(-2**n) * binomial(2**n, 2**(n-1)))
    else:
        c = c/2
        return 4 * big_precision(2**(-2**n) * binomial(2**n, 2**(n-1)+c))


def ddt_coeff_probability(m, n, c, precision=DEFAULT_HIGH_PRECISION):
    """Returns the probability that a coefficient of the DDT of a S-Box
    mapping m bits to n is equal to c.

    This probability is identical for random a permutation and a
    random non-bijective function.

    """
    big_precision = RealField(precision)
    if c % 2 == 1:
        return 0
    else:
        d = c/2
        if m >= 5 and n-m <= m/4:
            k = 2**(m-n-1)
            return RealNumber(exp(-k) * k**d / factorial(d))
        else:
            return RealNumber(binomial(2**(m-1), d) * 2**(-n*d) * (1 - 2**-n)**(2**(m-1)-d))


def bct_coeff_probability(m, n, c, precision=DEFAULT_HIGH_PRECISION):
    """Returns the probability that a coefficient of the BCT of an S-Box
    mapping m bits to n is equal to c.

    This probability is only defined for permutations. Thus, an error is raised if m != n.

    """
    big_precision = RealField(precision)
    if m != n:
        raise "the BCT is only defined when m==n"
    if c % 2 == 1:
        return RealNumber(0.0)
    B = big_precision(2**(n-1))
    A = big_precision(2**(2*n-2)-2**(n-1))
    p = big_precision(1/(2**n-1))
    q = p**2
    d = int(c/2)
    result = big_precision(0.0)
    base = big_precision(2**n-1)**(B + 2*A)
    for j1, j2 in itertools.product(range(0, d+1), range(0, d+1)):
        if 2*j1 + 4*j2 == c:
            # added = Integer(binomial(B, j1)) * RealNumber(p**j1) * RealNumber((1-p)**(B-j1)) * Integer(binomial(A, j2)) * RealNumber(q**(j2) * (1-q)**(A-j2))
            added = big_precision(binomial(B, j1)) * big_precision((2**n-2)**(B-j1)) * big_precision(binomial(A, j2)) * big_precision((2**(2*n)-2**(n+1))**(A-j2))
            if added > 0 and added < Infinity:
                result += added / base
    return result

        
def expected_max_ddt(m, n):
    result = RealNumber(0.0)
    cumul = [0]
    for k in range(1, 15):
        to_add = sum(ddt_coeff_probability(m, n, 2*z) for z in range(0, k+1))
        to_add = to_add**RealNumber((2**m-1) * (2**n-1))
        cumul.append(to_add)
    result = []
    for i in range(1, len(cumul)):
        result.append(cumul[i] - cumul[i-1])
    return result


def expected_max_lat(m, n):
    result = RealNumber(0.0)
    cumul = [0]
    for k in range(1, 50):
        to_add = sum(lat_coeff_probability_permutation(m, n, 2*z) for z in range(0, k+1))
        to_add = to_add**RealNumber((2**m-1) * (2**n-1))
        cumul.append(to_add)
    result = []
    for i in range(1, len(cumul)):
        result.append(cumul[i] - cumul[i-1])
    return result

def expected_max_lat_function(m, n):
    result = RealNumber(0.0)
    cumul = [0]
    for k in range(1, 70):
        to_add = sum(lat_coeff_probability_function(m, n, z) for z in range(0, k+1))
        to_add = to_add**RealNumber((2**m) * (2**n))
        cumul.append(to_add)
    result = []
    for i in range(1, len(cumul)):
        result.append(cumul[i] - cumul[i-1])
    return result
    
    
# !SUBSECTION! Aggregated information from the tables

def probability_of_max_and_occurrences(m, n, v_max, occurrences, proba_func, precision=DEFAULT_HIGH_PRECISION):
    """Returns the logarithm in base 2 of the probability that
    $(2^m-1)(2^n-1)$ trials of an experiment yielding output c with
    probability proba_func(m, n, c) will have a result equal to
    `v_max` at most `occurrences` times and be strictly smaller on all
    other trials.

    """
    try:
        big_precision = RealField(precision)
        p_strictly_smaller = sum(big_precision(proba_func(m, n, i, precision=precision)) for i in range(0, v_max))
        p_equal = big_precision(proba_func(m, n, v_max, precision=precision))
        result = big_precision(0)
        n_trials = (2**m-1) * (2**n-1)
        for occ in reversed(range(0, occurrences+1)):
            added = big_precision(binomial(n_trials, occ)) * (p_strictly_smaller)**((n_trials - occ)) * p_equal**(occ)
            result += added
        return result
    except:
        print("failure in sum: occ={}, result={}, added={}".format(occ, result, added))
        return 0


def anomaly_differential_uniformity(n, v_max):
    p_0 = RealNumber(ddt_coeff_probability(n, n, 0))
    p_non_zero = RealNumber(1-p_0)
    p_inf = RealNumber(sum(ddt_coeff_probability(n, n, c) for c in range(1, v_max+1)))
    proba_enough_zeroes_in_row = RealNumber(0.0)
    proba_bound_and_enough_zeroes = RealNumber(0.0)
    M = 2**(n-1)-1
    for i in range(2**(n-1)-1, 2**n-1):
        proba_enough_zeroes_in_row    += RealNumber(binomial(2**n-1, i)) * p_0**i * p_non_zero**(2**n-1-i)
        proba_bound_and_enough_zeroes += RealNumber(binomial(2**n-1, i)) * p_0**i * p_inf**(2**n-1-i)
    return (2**n-1)*float(log(proba_bound_and_enough_zeroes, 2) - log(proba_enough_zeroes_in_row, 2))


def anomaly_ddt(n, v_max, occ):
    p_0 = RealNumber(ddt_coeff_probability(n, n, 0))
    p_non_zero = RealNumber(1-p_0)
    p_inf = RealNumber(sum(ddt_coeff_probability(n, n, c) for c in range(1, v_max+1)))
    proba_enough_zeroes_in_row = RealNumber(0.0)
    proba_bound_and_enough_zeroes = RealNumber(0.0)
    M = 2**(n-1)-1
    for i in range(2**(n-1)-1, 2**n-1):
        proba_enough_zeroes_in_row    += RealNumber(binomial(2**n-1, i)) * p_0**i * p_non_zero**(2**n-1-i)
        proba_bound_and_enough_zeroes += RealNumber(binomial(2**n-1, i)) * p_0**i * p_inf**(2**n-1-i)
    return -(2**n-1)*float(log(proba_bound_and_enough_zeroes, 2) - log(proba_enough_zeroes_in_row, 2))
    

def table_anomaly(s, table, spec=None, precision=DEFAULT_HIGH_PRECISION):
    if table not in ["DDT", "LAT", "BCT"]:
        raise "The table-based anomaly is defined for the LAT, DDT and BCT. table={} is unknown.".format(table)
    else:
        proba_func = None
        n = int(log(len(s), 2))
        if table == "DDT":
            if spec == None:
                spec = differential_spectrum(s)
            proba_func = ddt_coeff_probability
        elif table == "LAT":
            if spec == None:
                spec = walsh_spectrum(s)
            proba_func = lat_coeff_probability_permutation
        else:
            if spec == None:
                spec = boomerang_spectrum(s)
            proba_func = bct_coeff_probability
        v_max = 0
        occurrences = 0
        for k in spec.keys():
            if abs(k) == v_max:
                occurrences += spec[k]
            elif abs(k) > v_max:
                v_max = abs(k)
                occurrences = spec[k]
        big_precision = RealField(precision)
        p = big_precision(probability_of_max_and_occurrences(n, n, v_max, occurrences, proba_func, precision=precision))
        return -p.log2()

    
def table_negative_anomaly(s, table, spec=None, precision=DEFAULT_HIGH_PRECISION):
    if table not in ["DDT", "LAT", "BCT"]:
        raise "The table-based negative anomaly is defined for the LAT, DDT and BCT. table={} is unknown.".format(table)
    else:
        proba_func = None
        n = int(log(len(s), 2))
        if table == "DDT":
            if spec == None:
                spec = differential_spectrum(s)
            proba_func = ddt_coeff_probability
        elif table == "LAT":
            if spec == None:
                spec = walsh_spectrum(s)
            proba_func = lat_coeff_probability_permutation
        else:
            if spec == None:
                spec = boomerang_spectrum(s)
            proba_func = bct_coeff_probability
        v_max = 0
        occurrences = 0
        for k in spec.keys():
            if abs(k) == v_max:
                occurrences += spec[k]
            elif abs(k) > v_max:
                v_max = abs(k)
                occurrences = spec[k]
        p_anomaly = probability_of_max_and_occurrences(n, n, v_max, occurrences, proba_func, precision=precision)
        p_equal = proba_func(n, n, v_max, precision=precision)
        p_strictly_smaller = sum(proba_func(n, n, i, precision=precision) for i in range(0, v_max))
        card = (2**n-1)**2
        big_precision = RealField(precision)
        p_precise_equal = big_precision(binomial(card, occurrences))*p_equal**occurrences*p_strictly_smaller**(card-occurrences)
        # p_precise_equal = min(p_precise_equal, big_precision(1.0))
        # p_anomaly       = min(p_anomaly, big_precision(1.0))
        return -big_precision(p_precise_equal-p_anomaly+1).log2()



# !SUBSECTION! Algebraic properies

@parallel
def algebraic_normal_form_coordinate(s, i):
    """Returns the algebraic normal form of the `i`th coordinate of the
    SBox s.

    """
    coordinate = BooleanFunction([(x >> i) & 1 for x in list(s)])
    return coordinate.algebraic_normal_form()


def algebraic_normal_form(s):
    """Returns the algebraic normal form of the coordinates of the S-Box
    `s`.

    If we let each of the n output bits of `s` be a Boolean function
    s_i so that $s(x) = s_{n-1}(x) s_{n-2}(x) ... s_{0}(x)$, then
    returns a list l such that l[i] is the ANF of s_i.

    """
    result = {}
    n = int(log(len(s), 2))
    outputs = algebraic_normal_form_coordinate([(s, b) for b in range(0, n)])
    for entry in outputs:
        result[entry[0][0][1]] = entry[1]
    return [result[k] for k in range(0, n)]



def algebraic_degree(s):
    """Returns the algebraic degree of `s`, i.e. the maximum algebraic
    degree of its coordinates.

    """
    anfs = algebraic_normal_form(s)
    result = 0
    for a in anfs:
        result = max(result, a.degree())
    return result


def degree_coordinates(s):
    """Returns a dictionnary `d` such that `d[k]==l` if and only if the
    function with LUT `s` has exactly `l` coordinates with algebraic degree
    `k`.

    """
    return [a.degree() for a in algebraic_normal_form(s)]
    

def degree_spectrum(s):
    """Returns a dictionnary `d` such that `d[k]==l` if and only if the
    function with LUT `s` has exactly `l` components with algebraic degree
    `k`.

    """
    anf = algebraic_normal_form(s)
    result = defaultdict(int)
    n = int(log(len(s), 2))
    for mask in range(1, 2**n):
            component = 0
            for i in range(0, n):
                if ((mask >> i) & 1) == 1:
                    component += anf[i]
            result[component.degree()] += 1
    return result


def hdim(s):
    """Returns the High Degree Indicator Matrix of `s`.

    Let n be such that `s` is the LUT of a function mapping F_2^n to
    itself. Then the HDIM is the n x n matrix with binary coefficients
    m_{i,j}, where

    m_{i,j} = \\sum_x (e_i \\cdot F(x)) (e_j \\cdot x)

    as defined in

    Perrin, Léo, and Aleksei Udovenko. "Algebraic insights into the
    secret Feistel network (full version)." Cryptology ePrint Archive
    (2016).

    """
    n = int(log(len(s), 2))
    H = [[0 for j in range(0, n)] for i in range(0, n)]
    b = [int(1 << i) for i in range(0, n)] # canonical basis
    for i, j in itertools.product(range(0, n), range(0, n)):
        H[i][j] = sum(scal_prod(b[i], s[x]) * scal_prod(b[j], x) for x in range(0, 2**n)) % 2
    return Matrix(GF(2), n, n, H)
    

    
# !SECTION! Inverting a LAT

def invert_lat(l):
    """Assuming that l is a valid function LAT, returns this function.

    Only works for functions with the same input and output length."""
    n = int(log(len(l), 2))
    return invert_lat_fast(l, n)




# !SECTION! Tests

if __name__ == '__main__':
    s = random_permutation(5)
    print(s)
    l = lat(s)
    s_prime = invert_lat(l)
    print(s_prime)

