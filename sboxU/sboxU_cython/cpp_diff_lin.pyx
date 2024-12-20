# -*-python-*- 
# Time-stamp: <2024-12-20 10:55:28>

import os
from libcpp cimport bool
from libcpp.vector cimport vector
from libc.stdint cimport int64_t, uint64_t
from sboxu_cpp cimport *
import os
from sage.all import Integer

BIG_SBOX_THRESHOLD = 128
DEFAULT_N_THREADS  = 2


cdef class PyFptFunction :
    cdef CppFptFunction fpt_function

    def __cinit__(self, int64_t p, int64_t t, int64_t u, vector[vector[int64_t]] indexes, vector[vector[int64_t]] values):
        self.fpt_function = CppFptFunction(p,t,u,indexes,values)

    def __call__(self, vector[int64_t] index):
        return self.fpt_function.eval(index)

    @property
    def p(self):
        return self.fpt_function.p
    @p.setter
    def p(self, p):
        self.fpt_function.p = p

    @property
    def t(self):
        return self.fpt_function.t
    @t.setter
    def t(self, t):
        self.fpt_function.t = t

    @property
    def u(self):
        return self.fpt_function.u
    @u.setter
    def u(self, u):
        self.fpt_function.u = u
    
    def fpt_ddt(self):
        return self.fpt_function.fpt_ddt()

    def fpt_differential_spectrum(self,n_threads=None):
        if n_threads == None:
            if self.p**self.t > BIG_SBOX_THRESHOLD:
                n_threads = DEFAULT_N_THREADS
            else:
                n_threads = 1
        return self.fpt_function.fpt_differential_spectrum_fast(n_threads)



def tab_vector_to_tab_int(tab_of_vectors, p):
    sol = [0] * len(tab_of_vectors)
    for i, vec in enumerate(tab_of_vectors):
        elt = 0
        for scalar in vec:
            elt = p*elt + scalar
        sol[i] = elt
    return sol

def tab_int_to_tab_vector(tab_of_ints, p, m):
    sol = [[0]*m for _ in range(len(tab_of_ints))]
    for i, int_rep in enumerate(tab_of_ints):
        temp = int_rep
        for j in range(m):
            sol[i][m - 1 - j] = temp%p
            temp //= p
    return sol    


# !SECTION! Differential properties

def differential_spectrum(s, p=None, n_threads=None):
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
    if len(s)%2 == 0:
        return differential_spectrum_fast(s, n_threads)
    if p == None:
        raise ValueError("Base prime p not given. Usage: differential_spectrum(f, p=None, n_threads=None)")        
    if not Integer(p).is_prime():
        raise ValueError("p is not a prime.")
    m = 1
    l = len(s)
    temp = p
    while temp < l:
        temp *= p
        m += 1
    if temp != l:
        raise ValueError(f"Table size is not a power of {p}.")
    if m == 1:
        return PyFptFunction(p, m, m, [[i] for i in range(p)], [[s[i]] for i in range(p)]).fpt_differential_spectrum(n_threads)
    if isinstance(s[0], int):
        s_linearized = tab_int_to_tab_vector(s, p, m)
        return PyFptFunction(p, m, m, tab_int_to_tab_vector(list(range(len(s))), p, m), s_linearized).fpt_differential_spectrum(n_threads)
    else:
        return PyFptFunction(p, m, m, tab_int_to_tab_vector(list(range(len(s))), p, m), s).fpt_differential_spectrum(n_threads)


def ddt(s, p=None):
    """Returns the DDT of the S-box `s`, i.e. a list of list `T` such that
    `T[a][b]==l` if and only if the equation s(x + a) + s(x) = b has
    exactly `l` solutions, where + denotes the addition in F_p^n (i.e. a
    bitwise XOR).
    If p is not given, it is assumed to be 2.
    """

    if len(s)%2 == 0:
        return ddt_cpp(s)
    if p == None:
        raise ValueError("Base prime p not given. Usage: ddt(f, p=None)")        
    if not Integer(p).is_prime():
        raise ValueError("p is not a prime.")
    m = 1
    l = len(s)
    temp = p
    while temp < l:
        temp *= p
        m += 1
    if temp != l:
        raise ValueError(f"Table size is not a power of {p}.")
    if m == 1:
        return PyFptFunction(p, m, m, [[i] for i in range(p)], [[s[i]] for i in range(p)]).fpt_ddt()
    if isinstance(s[0], int):
        s_linearized = tab_int_to_tab_vector(s, p, m)
        return PyFptFunction(p, m, m, tab_int_to_tab_vector(list(range(len(s))), p, m), s_linearized).fpt_ddt()
    else:
        return PyFptFunction(p, m, m, tab_int_to_tab_vector(list(range(len(s))), p, m), s).fpt_ddt()


def ddt_row(s, a):
    """Returns the row corresponding to input difference `a` in the DDT of
    the function with lookup table `s`.
    Only for p = 2.
    """
    return ddt_row_cpp(s, a)
    

def ortho_derivative(s):
    """Assuming that `s` corresponds to a quadratic APN function, returns
    the unique function `f` such that f(a) . (s(x+a)+s(x)+s(a)+s(0)) = 0
    for all a != 0 and for all x.

    If `s` is not a quadratic APN function, returns an empty list.

    """
    result = ortho_derivative_fast(s)
    if len(result) < 2:
        raise Exception("Not quadratic")
    else:
        return result
    

def is_ddt_row_max_smaller_than(s, a, u):
    """Returns true if and only if the maximum coefficient of row `a`
    of the DDT of the S-box `s` (i.e. the maximum value of DDT[a, b]
    for the given `a`) is at most equal to `u`.

    """
    return is_ddt_row_max_smaller_than_cpp(s, a, u)

    

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

def lat(s, p=None, n_threads=DEFAULT_N_THREADS):
    """Return the Linear Approximation Table of `s` (or its Walsh
    coefficients). It is a 2-dimensional array `T` such that

    T[a][b] = \sum_{x=0}^{p^n-1} (omega)^{ax + b s[x]}

    where `s` is the function whose LUT is given as the input to this
    function, omega is a complex pth root of unity (-1 when p == 2)
    and where n is the size of the input set in bits.
    If p is not given, it is assumed to be 2.

    """
    if len(s)%2 == 0 or p == 2:
        return lat_cpp(s)
    if p == None:
        raise ValueError("Base prime p not given and probably not 2... Usage: lat(f, p, n_threads)")
    if not Integer(p).is_prime():
        raise ValueError("p is not a prime.")
    m = 1
    l = len(s)
    temp = p
    while temp < l:
        temp *= p
        m += 1
    if temp != l:
        raise ValueError(f"Table size is not a power of {p}.")
    assert(p**m == l)
    return fpt_lat(s, p, m, n_threads)

def lat_column(s, b, p=None):
    """
    Return the column at index b of the LAT of `s`.
    If p is not given, it is assumed to be 2.
    """
    if len(s)%2 == 0 or p == 2:
        raise ValueError("Not implemented for sizes powers of two.")
    if p == None:
        raise ValueError("Base prime p not given and assumed != 2. Usage: lat_column(f, b, p=None)")
    if not Integer(p).is_prime():
        raise ValueError("p is not a prime.")
    m = 1
    l = len(s)
    temp = p
    while temp < l:
        temp *= p
        m += 1
    if temp != l:
        raise ValueError(f"Table size is not a power of {p}.")
    assert(p**m == l)
    return fpt_lat_column(s, p, m, b)
    
def lat_row(s, a, p=None):
    """
    Return the row at index a of the LAT of `s`.
    /!\ ONLY WORKS IF s IS A PERMUTATION
    If p is not given, it is assumed to be 2.
    """
    if len(s)%2 == 0 or p == 2:
        raise ValueError("Not implemented for sizes powers of two.")
    if p == None:
        raise ValueError("Base prime p not given and assumed != 2. Usage: lat_row(f, a, p=None)")
    if not Integer(p).is_prime():
        raise ValueError("p is not a prime.")
    m = 1
    l = len(s)
    temp = p
    while temp < l:
        temp *= p
        m += 1
    if temp != l:
        raise ValueError(f"Table size is not a power of {p}.")
    assert(p**m == l)
    if (len(set(s)) != len(s)):
        raise ValueError("lat_row is not implemented for functions that are not permutations. Try lat or lat_column.")
    return fpt_lat_row(s, p, m, a)

def lat_max(s, p=None, n_threads=DEFAULT_N_THREADS):
    """
    Return the linearity of `s`, i.e. the maximum of its coefficients.
    If p is not given, it is assumed to be 2.
    """
    if len(s)%2 == 0 or p == 2:
        raise ValueError("Not implemented for sizes powers of two.")
    if p == None:
        raise ValueError("Base prime p not given and assumed != 2. Usage: lat_max(f, p=None, n_threads=2)")
    if not Integer(p).is_prime():
        raise ValueError("p is not a prime.")
    m = 1
    l = len(s)
    temp = p
    while temp < l:
        temp *= p
        m += 1
    if temp != l:
        raise ValueError(f"Table size is not a power of {p}.")
    assert(p**m == l)
    return fpt_max_lat(s, p, m, n_threads)

def walsh_spectrum(s, p=None, epsilon=1e-5, n_threads=None):
    """
    Returns the walsh spectrum of s, a.k.a. a dict ws such that ws[value] is the number
    of entries of the LAT at value.
    If p is not given, it is assumed to be 2.
    If p = 2, the entries of the LAT are integers.
    If p != 2, they are float.
    """
    if n_threads == None:
        if len(s) > BIG_SBOX_THRESHOLD:
            n_threads = DEFAULT_N_THREADS
        else:
            n_threads = 1
    if len(s)%2 == 0 or p == 2:
        return walsh_spectrum_fast_cpp(s, n_threads)
    if p == None:
        raise ValueError("Base prime p not given and assumed != 2. Usage: walsh_spectrum(f, p=None, n_threads=None)")
    if not Integer(p).is_prime():
        raise ValueError("p is not a prime.")
    m = 1
    l = len(s)
    temp = p
    while temp < l:
        temp *= p
        m += 1
    if temp != l:
        raise ValueError(f"Table size is not a power of {p}.")
    assert(p**m == l)
    return fpt_walsh_spectrum(s, p, m, epsilon, n_threads)

def fourier_transform(l):
    return walsh_spectrum_coord(l)

def invert_lat_fast(t,n):
    return invert_lat_cpp(t, n)

def lat_zeroes_fast(l, n, n_threads=1):
    return lat_zeroes_cpp(l,n,n_threads)

def proj_lat_zeroes(l, n_threads=1):
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
