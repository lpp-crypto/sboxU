#!/usr/bin/sage

from sage.all import Matrix, GF, vector, log, randint

import itertools
import random
from utils import oplus
from sboxu_cpp import *
from linear import *
from display import *
from diff_lin import *
from hashlib import sha256
from collections import defaultdict

DEFAULT_N_THREADS  = 16



# !SECTION! Utils

def thickness(basis, N):
    """Returns the thickness of the vector space with basis `basis`, where
    this vector space is a subspace of $\F_2^{N+M}$ for the given `N` and
    for some $M$.

    The thickness is the rank of the projection of this space on its
    `N` bits of lowest weight.

    """
    MASK_N = sum(int(1 << i) for i in xrange(0, N))
    proj = [b & MASK_N for b in basis]
    return rank_of_vector_set(proj, 2*N)


def thickness_spectrum(s):
    """Returns a dictionary containing the thickness spectra of the
    function whose LUT is the list `s`.

    It first computes the Walsh zeroes of `s`, then look for vector
    spaces of dimension N in it. For each space, 

    """
    N = int(log(len(s), 2))
    z_s = lat_zeroes(s)
    minimal_bases = extract_bases(z_s, N, 2*N, number="fixed dimension")
    result = defaultdict(int)
    for basis in minimal_bases:
        result[thickness(basis, N)] += 1
    return dict(result)

def get_lat_zeroes_spaces(s, n_threads=DEFAULT_N_THREADS):
    """Returns a list containing the basis of each vector space of
    dimension n contained in the LAT zeroes of `s`.

    """
    return get_lat_zeroes_spaces_fast(s,
                                      int(log(len(s), 2)),
                                      int(n_threads))

    

# !SECTION! Linear/Affine Equivalence 


# !SUBSECTION! XOR equivalence 

def xor_equivalence(f, g):
    """Returns a pair [k0, k1] of integers such that, for all x:

    f[x] = g[x + k0] + k1,

    where "+" denotes the bitwise XOR.
    """
    N = int(log(len(f), 2))
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



# !SECTION! CCZ-equivalence to a permutation 

# !SUBSECTION! Extended Affine equivalence

def ea_equivalent_permutation_mappings(f):
    """Returns the list of all the linear functions L such that f(x) +
    L(x) is a permutation, where L is in fact the matrix
    representation of this function.

    """
    N = int(log(len(f), 2))
    mask = sum((1 << i) for i in xrange(0, N))
    z = lat_zeroes(f)
    spaces = extract_bases(z, N, 2*N, number="fixed dimension")
    result = []
    for b in spaces:
        proj_dict = {}
        v = linear_span(b)
        l = [-1 for x in xrange(0, 2**N)]
        for x in v:
            l[x & mask] = x >> N
        if -1 not in l:
            result.append(linear_function_lut_to_matrix(l).transpose())
    return result


# !SUBSECTION! General case

def ccz_equivalent_permutations(f, number="all permutations"):
    """Returns a list of permutations that are CCZ-equivalent to
    `f`. 

    The behaviour of the function depends on the value of `number`:

    - if it is set to "all permutations" then at least one permutation
      per affine-equivalence class of permutations is returned; but

    - if it is set to "just one" then the output contains at most one
      permutation, namely the first one found.

    """
    result = []
    N = int(log(len(f), 2))
    mask = sum(int(1 << i) for i in xrange(0, N))
    graph_f = [(x << N) | f[x] for x in xrange(0, 2**N)]
    bases = get_lat_zeroes_spaces(f)
    bases_by_dimensions = defaultdict(list)
    for b in bases:
        t1 = rank_of_vector_set([v >> N for v in b], N)
        t2 = rank_of_vector_set([v & mask for v in b], N)
        bases_by_dimensions[(t1 << N) | t2].append(b)
    for dim_pairs in itertools.product(bases_by_dimensions.keys(),
                                       bases_by_dimensions.keys()):
        t1, t2 = dim_pairs[0] >> N, dim_pairs[0] & mask
        u1, u2 = dim_pairs[1] >> N, dim_pairs[1] & mask
        if (t1 + u1) >= N and (t2 + u2) >= N:
            for b0, b1 in itertools.product(bases_by_dimensions[dim_pairs[0]],
                                            bases_by_dimensions[dim_pairs[1]]):
                if rank_of_vector_set(b0 + b1, 2*N) == 2*N:
                    L = Matrix(GF(2), 2*N, 2*N, [tobin(x, 2*N) for x in b0 + b1])
                    graph_g = [apply_bin_mat(word, L) for word in graph_f]
                    g = [-1 for x in xrange(0, 2**N)]
                    for word in graph_g:
                        x, y = word >> N, word & mask
                        g[x] = y
                    if -1 in g:
                        raise Exception("permutation ill defined!")
                    else:
                        result.append(g)
                        if number == "just one":
                            return result
    return result



# !SECTION! Exploring CCZ-class

def enumerate_ea_classes(f):
    """Returns a list containing at least one function from each of the
    EA-classes constituting the CCZ-class of `f`.

    Note that several functions in the same EA-class may be
    returned. Solving this problem is in fact an open *research*
    problem.

    """
    N = int(log(len(f), 2))
    mask = sum(int(1 << i) for i in xrange(0, N))
    graph_f = [(x << N) | f[x] for x in xrange(0, 2**N)]
    bases = get_lat_zeroes_spaces(f)
    result = []
    for b in bases:
        L_map = FastLinearMapping(get_generating_matrix(b, 2*N).transpose())
        graph_g = [L_map(word) for word in graph_f]
        g = [-1 for x in xrange(0, 2**N)]
        for word in graph_g:
            x, y = word >> N, word & mask
            g[x] = y
        if -1 in g:
            raise Exception("permutation ill defined!")
        else:
            result.append(g)
    return result
                        


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



# !SUBSECTION!  Test CCZ-equivalent permutation

def test_ea_permutations():
    for N in [4, 5]:
        F = GF(2**N, name="a")
        inv = [(F.fetch_int(x)**(2**N-2)).integer_representation()
               for x in xrange(0, 2**N)]
        print "== ", N
        for L in ea_equivalent_permutation_mappings(inv):
            print L.str(), "\n"

def test_ccz_permutations(number="all permutations"):
    N = 6
    F = GF(2**N, name="a")
    # generating the Kim mapping
    kim = []
    for x_i in xrange(0, 2**N):
        x = F.fetch_int(x_i)
        y = x**3 + x**10 + F.gen()*x**24
        kim.append(y.integer_representation())
    permutations = ccz_equivalent_permutations(kim, number=number)
    for i, p in enumerate(permutations):
        print "{:2d} {} {} {}".format(
            i,
            is_permutation(p),
            pretty_spectrum(differential_spectrum(p)),
            pretty_vector(p)            
        )
    print "total: {}".format(len(permutations))

def test_enumerate_ea():
    N = 6
    F = GF(2**N, name="a")
    # generating the Kim mapping
    kim = []
    for x_i in xrange(0, 2**N):
        x = F.fetch_int(x_i)
        y = x**3 + x**10 + F.gen()*x**24
        kim.append(y.integer_representation())
    classes = enumerate_ea_classes(kim)
    for f in classes:
        print algebraic_degree(f), pretty_spectrum(thickness_spectrum(f))
    print "total: ", len(classes)
    

# !SECTION! Running tests

if __name__ == '__main__':
    # test_ea_permutations()
    # test_ccz_permutations(number="just one")
    test_enumerate_ea()
    
    # import sys
    # N = int(sys.argv[1])
    # print "=== Linear Equivalence ==="
    # test_le_equivalence(N, verbose=True)
    # print "\n=== Linear Representative ==="
    # test_le_repr(N, verbose=False)
    # print "\n=== Affine Equivalence ==="
    # test_ae_equivalence(N, verbose=True)
