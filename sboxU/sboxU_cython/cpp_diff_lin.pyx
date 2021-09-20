# -*-python-*- 
# Time-stamp: <2021-09-20 16:25:20 lperrin>

from sboxu_cpp cimport *

BIG_SBOX_THRESHOLD = 128
DEFAULT_N_THREADS  = 2


# !SECTION! Differential properties

def differential_spectrum(s, n_threads=None):
    """Returns the differential spectrum of the S-box `s`, i.e. a
    dictionnary `d` such that `d[k]==l` if there are `l` occurrences of
    the value `k` in the DDT of `s`.

    The optional parameter `n_threads` can be used to specify the
    number of threads used to perform this computation.

    """
    if n_threads == None:
        if len(s) > BIG_SBOX_THRESHOLD:
            n_threads = DEFAULT_N_THREADS
        else:
            n_threads = 1
    return differential_spectrum_fast(s, n_threads)


def ddt(s):
    """Returns the DDT of the S-box `s`, i.e. a list of list `T` such that
    `T[a][b]==l` if and only if the equation s(x + a) + s(x) = b has
    exactly `l` solutions, where + denotes the addition in F_2^n (i.e. a
    bitwise XOR).

    """
    return ddt_cpp(s)


def ortho_derivative(s):
    """Assuming that `s` corresponds to a quadratic APN function, returns
    the unique function `f` such that f(a) . (s(x+a)+s(x)+s(a)+s(0)) = 0
    for all a != 0 and for all x.

    If `s` is not a quadratic APN function, returns an empty list.

    """
    return ortho_derivative_fast(s)
    

def is_differential_uniformity_smaller_than(s, u):
    """Returns true if and only if the differential uniformity of the
    S-box `s` (i.e. the maximum value of DDT[a, b] for a > 0) is at
    most equal to `u`.

    """
    return is_differential_uniformity_smaller_than_cpp(s, u)



# !SUBSECTION! c-differential properties

def log_exp_table(F):
    """Returns a pair of tables of integers. Let `a` be a root of the
    polynomial generator of the multiplicative

    - The first, L, is a logarithm table. It is such that a**L[x] = x,
      and L[0] = -1.
    - The second, E, is the inverse of L: E[L[x]] = x.

    """
    a = F.gen()
    x = a
    N = F.degree()
    E = [1]
    for i in range(1, 2**N-1):
        E.append(x.integer_representation())
        x *= a
    L = [2**N]
    for i in range(1, 2**N):
        L.append(E.index(i))
    return [L, E]


def c_differential_spectra(s, F, l_table=None, e_table=None):
    """Returns a vector v such that v[c] is a map of the form {k : n_k}
    such that s[x^a] ^ c s[x] = b has exactly k solutions x for n_k
    different pairs (a,b).

    """
    if l_table == None or e_table == None:
        l_table, e_table = log_exp_table(F)
    return c_differential_spectra_cpp(s, l_table, e_table)



# !SECTION! Linear properties (Walsh/Fourier)

def lat(s):
    """Return the Linear Approximation Table of `l` (or its Walsh
    coefficients). It is a 2-dimensional array `T` such that

    T[a][b] = \sum_{x=0}^{2^n-1} (-1)^{ax + b s[x]}

    where `s` is the function whose LUT is given as the input to this
    function, and where n is the size of the input set in bits.

    """
    return lat_cpp(s)
    
def walsh_spectrum(l, n_threads=None):
    if n_threads == None:
        if len(l) > BIG_SBOX_THRESHOLD:
            n_threads = DEFAULT_N_THREADS
        else:
            n_threads = 1
    return walsh_spectrum_fast_cpp(l, n_threads)

def fourier_transform(l):
    return walsh_spectrum_coord(l)

def invert_lat_fast(t,n):
    return invert_lat_cpp(t, n)

def lat_zeroes_fast(l, n, n_threads=None):
    return lat_zeroes_cpp(l,n,n_threads)

def proj_lat_zeroes(l, n_threads=None):
    if n_threads == None:
        if len(l) > BIG_SBOX_THRESHOLD:
            n_threads = DEFAULT_N_THREADS
        else:
            n_threads = 1
    return projected_lat_zeroes_cpp(l,n_threads)


# !SECTION! Boomerang properties

def bct(l):
    return bct_cpp(l)

def boomerang_spectrum(l, n_threads=None):
    if n_threads == None:
        if len(l) > BIG_SBOX_THRESHOLD:
            n_threads = DEFAULT_N_THREADS
        else:
            n_threads = 1
    return bct_spectrum_fast_cpp(l, n_threads)
