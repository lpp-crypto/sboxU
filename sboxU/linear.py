#!/usr/bin/sage

from sage.all import Matrix, GF, vector, log, randint

import itertools
import random
from sboxu_cpp import *
from display import pretty_spectrum
from hashlib import sha256
from collections import defaultdict


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
    n = mat.ncols()
    x = Matrix(GF(2), n, 1, tobin(x, n))
    y = mat * x
    return frombin(y.T[0][:mat.nrows()])


def apply_bin_mat_lsb_first(x, mat):
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



# !SUBSECTION! Interacting with vector spaces and their supersets

def rank_of_vector_set(V, n):
    """Returns the rank of the matrix obtained by "stacking" the n-bit
    binary representation of the numbers in V.

    """
    M = Matrix(GF(2), len(V), n, [tobin(x, n) for x in V])
    return M.rank()



def extract_basis(v, N):
    """Returns a subset of v such that v is included that the span of
    these elements is at least as big as v. In particular, if v is a
    vector space, it returns a base.

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
    contains basis.

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


def span(basis, with_zero=True):
    result = []
    if with_zero:
        result.append(0)
    visited = defaultdict(int)
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


# !SECTION! Linear/Affine Equivalence 


# !SUBSECTION! XOR equivalence 

def xor_equivalence(f, g):
    """Returns a pair [k0, k1] of integers such that, for all x:

    f[x] = g[x + k0] + k1,

    where "+" denotes the bitwise XOR.
    """
    for k0 in xrange(0, 2**N):
        k1 = oplus(f[0], g[k0])
        good = True
        for x in xrange(1, 2**N):
            if oplus(f[x], g[oplus(k0, x)]) != k1:
                good = False
                break
        if good:
            return [k0, k1]
    return []


# !SUBSECTION! Linear equivalence

def linear_equivalence(f, g):
    """Returns, if it exists, the pair A, B of matrix such that, for all x:

    f(x) = (B o g o A)(x),

    where "o" denotes functional composition. If no such linear
    permutations exist, returns an empty list.

    Internally calls a function written in C++ for speed which
    implements the "Linear Equivalence (LE)" algorithm from

    Alex Biryukov, Christophe De Canniere, An Braeken, and Bart
    Preneel (2003).  "A Toolbox for Cryptanalysis: Linear and Affine
    Equivalence Algorithms", Advances in Cryptology -- EUROCRYPT 2003,
    Lecture Notes in Computer Science 2656, E. Biham (ed.),
    Springer-Verlag, pp. 33--50, 2003.

    """
    if len(f) != len(g):
        raise "f and g are of different dimensions!"
    if (f[0] == 0 and g[0] != 0) or (f[0] != 0 and g[0] == 0):
        return []
    result = linear_equivalence_fast(f, g)
    if len(result) == 0:
        return result
    A = linear_function_lut_to_matrix(result[0])
    B = linear_function_lut_to_matrix(result[1])
    return A, B


# !SUBSECTION! Affine equivalence 

def hash_sbox(f):
    """Returns a 64-char string obtained by hashing the base 10
    representation of each entry of the lookup table f with SHA-256.

    """
    hf = sha256()
    for x in f:
        hf.update(str(x))
    return hf.hexdigest()
    

def affine_equivalence(f, g):
    """Returns, if it exists, the tuple A, a, B, b where A and B are
    matrices and where a and b are integers such that, for all x:

    f(x) = (B o g o A)(x + a) + b,

    where "o" denotes functional composition and "+" denotes XOR. If
    no such affine permutations exist, returns an empty list.

    Internally calls a function written in C++ for speed which returns
    the "Linear Representative" using an algorithm from

    Alex Biryukov, Christophe De Canniere, An Braeken, and Bart
    Preneel (2003).  "A Toolbox for Cryptanalysis: Linear and Affine
    Equivalence Algorithms", Advances in Cryptology -- EUROCRYPT 2003,
    Lecture Notes in Computer Science 2656, E. Biham (ed.),
    Springer-Verlag, pp. 33--50, 2003.

    """
    if len(f) != len(g):
        raise "f and g are of different dimensions!"
    table_f = defaultdict(list)
    table_c = defaultdict(int)
    for c in xrange(0, len(f)):
        f_c = le_class_representative([oplus(f[x], c) for x in xrange(0, len(f))])
        d = hash_sbox(f_c)
        table_f[d] = f_c
        table_c[d] = c
    rs = []
    a = -1
    b = -1    
    for c in xrange(0, len(f)):
        g_c = le_class_representative([g[oplus(x, c)] for x in xrange(0, len(f))])
        d = hash_sbox(g_c)
        if d in table_c.keys():
            a = c
            b = table_c[d]
            rs = g_c
            break
    if a == -1:
        return []
    l_f = linear_equivalence([oplus(f[x], b) for x in xrange(0, len(f))],
                             rs)
    A_f, B_f = l_f[0], l_f[1]
    l_g = linear_equivalence([g[oplus(x, a)] for x in xrange(0, len(f))],
                             rs)
    A_g, B_g = l_g[0], l_g[1]
    A = A_g.inverse() * A_f
    B = B_f * B_g.inverse()
    a = apply_bin_mat(a, A.inverse())
    return [A, a, B, b]


# !SECTION! Tests

def print_result(n_valid, n_tested):
    verdict = "[success]"
    if (n_valid != n_tested):
        verdict = "[FAIL]"
    print "{} success rate: {}/{} = {:.03f}".format(
        verdict,
        n_valid,
        n_tested,
        float(n_valid)/n_tested)

# !SUBSECTION! Linear equivalence 

def check_linear_equivalence(f, g, A, B):
    """Returns True if and only if f = B o g o A."""
    for x in xrange(0, 2**N):
        y = apply_bin_mat(g[apply_bin_mat(x, A)], B)
        if y != f[x]:
            return False
    return True


def test_le_equivalence(N, verbose=False):
    from timeit import default_timer
    false_negatives = 0
    false_positives = 0
    n_tested = 500
    for i in xrange(0, n_tested):
        # checking if linearly equivalent permutations are identified
        g = random_permutation(N)
        A = rand_linear_permutation(N)
        B = rand_linear_permutation(N)
        f = [apply_bin_mat(g[apply_bin_mat(x, A)], B) for x in xrange(0, 2**N)]
        computation_start = default_timer()
        result = linear_equivalence(f, g)
        computation_end = default_timer()
        elapsed = computation_end - computation_start
        if len(result) > 1:
            if A != result[0] or B != result[1]:
                # if f is self-linear equivalent, matrices other than
                # (A,B) can be correct. We check for those.
                if check_linear_equivalence(f, g, result[0], result[1]):
                    if verbose:
                        print "[success]     LE {:0.4f}s (other matrices found)".format(elapsed)
                else:
                    false_negatives += 1
                    if verbose:
                        print "[FAIL]     LE {:0.4f}s (wrong matrices found)".format(elapsed)
            else:
                if verbose:
                    print "[success]     LE {:0.4f}s".format(elapsed)
        else:
            false_negatives += 1
            if verbose:
                print "[FAIL]     LE {:0.4f} (nothing found)".format(
                    computation_end - computation_start)
        # checking if non-linearly equivalent functions are identified
        g = random_permutation(N)
        result = linear_equivalence(f, g)
        if len(result) == 0:
            if verbose:
                print "[success] non-LE {:0.4f}".format(elapsed)
        else:
            if check_linear_equivalence(f, g, result[0], result[1]):
                if verbose:
                    print "[success] act.LE {:0.4f}".format(elapsed)
            else:
                false_positives += 1
                if verbose:
                    "[FAIL] matrices found for non-LE permutations"
    print "* testing if LE functions are identified correctly (with correct linear permutations)"
    print_result(n_tested-false_negatives, n_tested)
    print "* testing if NON-LE functions are identified correctly"
    print_result(n_tested-false_positives, n_tested)



# !SUBSECTION! LE representative 

def test_le_repr(N, verbose=False):
    import diff_lin
    from timeit import default_timer
    n_valid = 0
    n_tested = 50
    print "* Testing whether f and le_class_representative(f) are LE"
    for i in xrange(0, n_tested):
        f = random_permutation(N)
        computation_start = default_timer()
        g = le_class_representative(f)
        computation_end = default_timer()
        if len(linear_equivalence(f, g)) != 0:
            n_valid += 1
            if verbose:
                print "[success] {:0.4f}s".format(computation_end - computation_start)
        else:
            print "[FAIL]"
    print_result(n_valid, n_tested)
    print "* testing whether two linear equivalent functions have the same representative"
    n_valid = 0
    for i in xrange(0, n_tested):
        f = random_permutation(N)
        A = rand_linear_permutation(N)
        B = rand_linear_permutation(N)
        g = [apply_bin_mat(f[apply_bin_mat(x, A)], B) for x in xrange(0, 2**N)]
        rs_f = le_class_representative(f)
        rs_g = le_class_representative(g)
        identical = True
        for x in xrange(0, 2**N):
            if rs_f[x] != rs_g[x]:
                identical = False
                break
        if identical:
            n_valid += 1
            if verbose:
                print "[success]"
        else:
            if verbose:
                print "[FAIL] representatives don't match"
                print rs_f, pretty_spectrum(diff_lin.differential_spectrum(rs_f))
                print rs_g, pretty_spectrum(diff_lin.differential_spectrum(rs_g))
    print_result(n_valid, n_tested)
    

# !SUBSECTION! Affine equivalence
    
def check_affine_equivalence(f, g, A, a, B, b):
    """Checks whether f(x) = (B o g o A)(x + a) + b"""
    for x in xrange(0, 2**N):
        y = oplus(x, a)
        y = apply_bin_mat(y, A)
        y = g[y]
        y = apply_bin_mat(y, B)
        y = oplus(y, b)
        if y != f[x]:
            return False
    return True


def test_ae_equivalence(N, verbose=False):
    from timeit import default_timer
    false_negatives = 0
    false_positives = 0
    n_tested = 10
    for i in xrange(0, n_tested):
        # checking if linearly equivalent permutations are identified
        f = random_permutation(N)
        A = rand_linear_permutation(N)
        a = randint(0, 2**N-1)
        B = rand_linear_permutation(N)
        b = randint(0, 2**N-1)
        g = [oplus(apply_bin_mat(f[apply_bin_mat(oplus(x, a), A)], B), b)
             for x in xrange(0, 2**N)]
        computation_start = default_timer()
        result = affine_equivalence(f, g)
        computation_end = default_timer()
        elapsed = computation_end - computation_start
        if len(result) > 1:
            if not check_affine_equivalence(f, g, result[0], result[1], result[2], result[3]):
                false_negatives += 1
                if verbose:
                    print "[FAIL] wrong affine permutations"
            else:
                if verbose:
                    print "[success]     AE {:0.4f}".format(elapsed)

        else:
            false_negatives += 1
            if verbose:
                print "[FAIL]     AE {:0.4f}s (nothing found)".format(elapsed)
        # checking if non-affine equivalent functions are identified
        g = random_permutation(N)
        result = affine_equivalence(f, g)
        if len(result) == 0:
            if verbose:
                print "[success] non-AE {:0.4f}s".format(elapsed)
        else:
            if check_affine_equivalence(f, g, result[0], result[1], result[2], result[3]):
                if verbose:
                    print "[success] act.AE {:0.4f}".format(elapsed)
            else:
                false_positives += 1
                if verbose:
                    "[FAIL] matrices found for non-LE permutations"
    print "* testing if AE functions are identified correctly (with correct affine permutations)"
    print_result(n_tested-false_negatives, n_tested)
    print "* testing if NON-LE functions are identified correctly"
    print_result(n_tested-false_positives, n_tested)

    

# !SECTION! Running tests

if __name__ == '__main__':
    import sys
    N = int(sys.argv[1])
    # print "=== Linear Equivalence ==="
    # test_le_equivalence(N, verbose=True)
    # print "\n=== Linear Representative ==="
    # test_le_repr(N, verbose=False)
    print "\n=== Affine Equivalence ==="
    test_ae_equivalence(N, verbose=True)
