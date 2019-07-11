#!/usr/bin/sage

from sage.all import Matrix, GF, vector, log, randint

import itertools
import random
from utils import oplus
from sboxu_cpp import *
from display import pretty_spectrum
from diff_lin import lat_zeroes
from hashlib import sha256
from collections import defaultdict

DEFAULT_N_THREADS  = 16


# !SECTION! Utils for dealing with linear functions

def tobin(x, n):
    return [(x >> i) & 1 for i in reversed(xrange(0, n))]

def frombin(v):
    y = 0
    for i in xrange(0, len(v)):
        y = (y << 1) | int(v[i])
    return y


# !SUBSECTION! Generating random functions

def rand_linear_permutation(n, density=0.5):
    """Returns the matrix representation of a linear permutation of
    {0,1}^n.

    This is done by generating random n x n binary matrices until one
    with full rank is found. As a random binary matrix has full rank
    with probability more than 1/2, this method is fine.

    """
    while True:
        result = [[0 for j in xrange(0, n)] for i in xrange(0, n)]
        for i, j in itertools.product(xrange(0, n), repeat=2):
            if random.random() < density:
                result[i][j] = 1
        result = Matrix(GF(2), n, n, result)
        if result.rank() == n:
            return result                
        
def rand_linear_function(m, n, density=0.5):
    """Returns the matrix representation of a linear function mapping
    {0,1}^m to {0,1}^m.

    """
    result = [[0 for j in xrange(0, m)] for i in xrange(0, n)]
    for i, j in itertools.product(xrange(0, n), xrange(0, m)):
        if random.random() < density:
            result[i][j] = 1
    return Matrix(GF(2), n, m, result)



# !SUBSECTION! Function composition

def apply_bin_mat(x, mat):
    """Interprets `x` as a binary vector corresponding to the binary
    representation of its value and multiplies it by the binary matrix
    `mat`.

    The result is returned as an integer whose binary representation
    is said product.

    """
    n = mat.ncols()
    x = Matrix(GF(2), n, 1, tobin(x, n))
    y = mat * x
    return frombin(y.T[0][:mat.nrows()])


def apply_bin_mat_lsb_first(x, mat):
    """Same as `apply_bin_mat` except that the order of the bits in `x`
    and in the output are reversed.
    
    """
    n = mat.nrows()
    bin_x = tobin(x, n)
    bin_x.reverse()
    x = Matrix(GF(2), n, 1, bin_x)
    y = mat * x
    bin_y = list(y.T[0][:n])
    bin_y.reverse()
    return frombin(bin_y)


def apply_bit_permutation(x, p):
    """Let $x = \sum x_i 2^i$ and p = [p_0, ... p_{n-1}]. Returns
    $y = \sum x_{p_i} 2^i$.

    """
    n = len(p)
    bin_x = [(x >> i) & 1 for i in xrange(0, n)]
    return sum(int(bin_x[p[i]] << i) for i in xrange(0, n))

def swap_halves(x, n):
    """If n=2k, swaps the first k bits of x with its last k bits."""
    if n % 2 == 1:
        raise Exception("Can't cut a word of {} bits in 2!".format(n))
    mask = sum(1 << i for i in xrange(0, n/2))
    l = x >> (n/2)
    r = x & mask
    return (r << (n/2)) | l


# !SUBSECTION! Linear functions and their LUT

def linear_function_lut_to_matrix(l):
    """Turns the look up table of a linear function into the
    corresponding binary matrix."""
    n = int(log(len(l), 2))
    result = []
    for i in xrange(0, n):
        line = [(int(l[1 << (n-1-j)]) >> (n-1-i)) & 1 for j in xrange(0, n)]
        result.append(line)
    return Matrix(GF(2), n, n, result)


def linear_function_matrix_to_lut(mat):
    """Returns the codebook of the matrix mat."""
    result = [0 for x in xrange(0, 2**mat.ncols())]
    for x in xrange(0, 2**mat.ncols()):
        x_vec = tobin(x, mat.ncols())
        y = frombin(mat * vector(x_vec))
        result[x] = y
    return result


def partial_linear_permutation_to_full(v, n):
    """Returns a matrix corresponding to a linear permutation of n bits
    mapping any x < len(v) to v[x].

    """
    if len(v) > 2**n:
        raise "n should be at least as big as the dimension of the state"
    if v[0] != 0:
        raise "v should start with 0"
    u = int(log(len(v), 2))
    if len(v) != 2**u:
        raise "size of v should be a power of 2"
    basis = []
    rebuilt_space = [0]
    # finding a basis of v
    for new_vector in v:
        if new_vector not in rebuilt_space:
            basis.append(new_vector)
            new_rebuilt_space = list(rebuilt_space)
            for x in rebuilt_space:
                new_rebuilt_space.append(oplus(x, new_vector))
            rebuilt_space = new_rebuilt_space
        if len(basis) == u:
            break
    # completing the basis
    while len(rebuilt_space) < 2**n:
        new_vector = 0
        while new_vector in rebuilt_space:
            new_vector += 1
        basis.append(new_vector)
        new_rebuilt_space = list(rebuilt_space)
        for x in rebuilt_space:
            new_rebuilt_space.append(oplus(x, new_vector))
        rebuilt_space = new_rebuilt_space
    result = Matrix(GF(2), n, n, [
        [(b >> (n - 1 - i)) & 1 for i in xrange(0, n)]
        for b in reversed(basis)]).transpose()
    check = linear_function_matrix_to_lut(result)
    if check[0:len(v)] == v:
        return result
    else:
        raise "no such matrix"



# !SUBSECTION! Vector/affine space extraction


def extract_bases(z,
                  dimension,
                  word_length,
                  n_threads=DEFAULT_N_THREADS,
                  number="all dimensions"):
    """Returns a list containing the Gaussian Jacobi basis of each vector
    space of dimension `dimension` that is contained in the list `z` of
    integers intepreted as elements of $\F_2^n$ where $n$ is equal to
    `word_length`.

    The number of threads to use can be specified using the argument
    `n_threads`.

    It can have 3 different behaviours depending on the value of the
    argument `number`:

    - if it is "just one" then it will stop as soon as it has found a
      vector space of the desired dimension and return its basis;

    - if it is "fixed dimension" then it will return all vector spaces
      with exactly the given dimension, even if they are subspaces of
      a larger space contained in `z`; and

    - if it is "all dimensions" then it will return all vector spaces
      of dimensions at least `dimension`. If a larger vector space is
      found, its bases will be return and its subspaces will be
      ignored.

    """
    if number not in ["all dimensions", "fixed dimension", "just one"]:
        raise Exception("Unknown value for parameter `number` in extract_bases:" + number)
    result = extract_bases_fast(z,
                                int(dimension),
                                int(word_length),
                                int(n_threads),
                                str(number))
    if number == "all dimensions":
        # in the case where we have obtained larger spaces, we remove
        # their subspaces from the list
        bigger_spaces = [b for b in result if len(b) > dimension]
        if len(bigger_spaces) == 0:
            # nothing to do as there is no bigger space
            return result
        else:
            new_result = list(bigger_spaces)
            for b in result:
                if len(b) == dimension:
                    is_included_in_bigger = False
                    for big_b in bigger_spaces:
                        is_included = True
                        for v in b:
                            if v not in big_b:
                                is_included = False
                                break
                        if is_included:
                            is_included_in_bigger = True
                            break
                    if not is_included_in_bigger:
                        new_result.append(b)
            return new_result
    else:
        return result
                            
            

def extract_affine_bases(z,
                         dimension,
                         word_length,
                         n_threads=DEFAULT_N_THREADS,
                         number="all dimensions"):
    """Returns a list containing the Gaussian Jacobi basis of each affine
    space of dimension `dimension` that is contained in the list `z` of
    integers intepreted as elements of $\F_2^n$ where $n$ is equal to
    `word_length`.

    The number of threads to use can be specified using the argument
    `n_threads`.

    It can have 3 different behaviours depending on the value of the
    argument `number`:

    - if it is "just one" then it will stop as soon as it has found an
      affine space of the desired dimension and return its basis;

    - if it is "fixed dimension" then it will return all affine spaces
      with exactly the given dimension, even if they are affine
      subspaces of a larger space contained in `z`; and

    - if it is "all dimensions" then it will return all vector spaces
      of dimensions at least `dimension`. If a larger vector space is
      found, its bases will be return and its subspaces will be
      ignored.

    """
    if number not in ["all dimensions", "fixed dimension", "just one"]:
        raise Exception("Unknown value for parameter `number` in extract_affine_bases:" + number)
    result = extract_affine_bases_fast(z,
                                     int(dimension),
                                     int(word_length),
                                     int(n_threads),
                                     str(number))
    if number == "all dimensions":
        # in the case where we have obtained larger spaces, we remove
        # their subspaces from the list
        bigger_affine = [[oplus(b[0], x) for x in linear_span(b[1:])]
                         for b in result if len(b) > dimension + 1]
        if len(bigger_affine) == 0:
            # nothing to do as there is no bigger space
            return result
        else:
            new_result = list([b for b in result if len(b) > dimension+1])
            for b in result:
                if len(b) == dimension+1:
                    aff = [oplus(b[0], v) for v in linear_span(b[1:])]
                    is_included_in_bigger = False
                    for big_space in bigger_affine:
                        is_included = True
                        for v in aff:
                            if v not in big_space:
                                is_included = False
                                break
                        if is_included:
                            is_included_in_bigger = True
                            break
                    if not is_included_in_bigger:
                        new_result.append(b)
            return new_result
    else:
        return result
    


# !SUBSECTION!  Vector space bases and their properties


def rank_of_vector_set(V, n):
    """Returns the rank of the matrix obtained by "stacking" the n-bit
    binary representation of the numbers in V.

    """
    M = Matrix(GF(2), len(V), n, [tobin(x, n) for x in V])
    return M.rank()



def extract_basis(v, N):
    """Returns a subset of v such that the span of these elements is at
    least as big as v. In particular, if v is a vector space, it returns a
    base.

    """
    dim = rank_of_vector_set(v, N)
    i = 0
    basis = []
    while i < len(v) and v[i] == 0:
        i += 1
    if i == len(v):
        return []
    basis.append(v[i])
    if dim == 1:
        return basis
    i += 1
    r = Matrix(GF(2), 1, N, [tobin(x, N) for x in basis]).rank()
    while r < dim and i < len(v):
        new_basis = basis + [v[i]]
        new_r = Matrix(GF(2), len(new_basis), N, [tobin(x, N) for x in new_basis]).rank()
        if new_r == dim:
            return new_basis
        elif new_r > r:
            basis = new_basis
            r = new_r
        i += 1
    return []


def complete_basis(basis, N):
    """Returns a list of length N spanning the space 0..2^N-1 which
    contains the list of integers `basis`. Assumes that the elements
    of `basis` are linearly independent.

    """
    r = len(basis)
    for i in xrange(1, 2**N):
        new_basis = basis + [i]
        new_r = Matrix(GF(2), len(new_basis), N, [tobin(x, N) for x in new_basis]).rank()
        if new_r == N:
            return new_basis
        elif new_r > r:
            basis = new_basis
            r = new_r
    return []


def get_generating_matrix(basis, N):
    """Returns an NxN binary matrix M such that M*(1 << i) = basis[i] for
    all i < len(basis) and such that M has full rank.

    """
    b = complete_basis(basis, N)
    return Matrix(GF(2), N, N, [
        [(b[i] >> j) & 1 for j in reversed(xrange(0, N))]
        for i in xrange(0, N)
    ]).transpose()


def linear_span(basis, with_zero=True):
    result = []
    if with_zero:
        result.append(0)
    visited = defaultdict(int)
    visited[0] = 1
    for i in xrange(1, 2**len(basis)):
        x = 0
        for j in xrange(0, len(basis)):
            if (i >> j) & 1 == 1:
                x = oplus(x, basis[j])
        if visited[x] != 1:
            result.append(x)
            visited[x] = 1
    return result

    

# !SUBSECTION! Easy interaction with finite fields

def mult_ff(x, y, F):
    return (F.fetch_int(x) * F.fetch_int(y)).integer_representation()


def div_ff(x, y, F):
    return (F.fetch_int(x) / F.fetch_int(y)).integer_representation()


def pow_ff(x, a, F):
    return (F.fetch_int(x)**a).integer_representation()


