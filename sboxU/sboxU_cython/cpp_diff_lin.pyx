# -*-python-*- 
# Time-stamp: <2021-08-11 15:43:20 lperrin>

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

# !SECTION! Linear properties (Walsh/Fourier)

def lat(l):
    return lat_cpp(l)
    
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
