#!/usr/bin/sage

from sage.all import Matrix, GF, vector, log, randint

import itertools
import random
from hashlib import sha256
from collections import defaultdict

from .utils import oplus
from .sboxU_cython import *
from .linear import *
from .display import *
from .diff_lin import *


DEFAULT_N_THREADS  = 16



# !SECTION! Utils


# !SUBSECTION! Basic CCZ-class invariant


def gamma_rank(f):
    """Returns the Gamma-rank of the function with LUT f.

    The Gamma-rank is the rank of the 2^{2n} \\times 2^{2n} binary
    matrix M defined by

    M[x][y] = 1 if and only if x + y \\in \\Gamma,

    where \\Gamma is the codebook of f, i.e.

    \\Gamma = \\{ (x, f(x)), x \\in \\F_2^n  \\} ~.

    """
    n = int(log(len(f), 2))
    dim = 2**(2*n)
    gamma = [(x << n) | f[x] for x in range(0, 2**n)]
    mat_content = []
    for x in range(0, dim):
        row = [0 for j in range(0, dim)]
        for y in gamma:
            row[oplus(x, y)] = 1
        mat_content.append(row)
    mat_gf2 = Matrix(GF(2), dim, dim, mat_content)
    return mat_gf2.rank()


def delta_rank(f):
    """Returns the Gamma-rank of the function with LUT f.

    The Gamma-rank is the rank of the 2^{2n} \\times 2^{2n} binary
    matrix M defined by

    M[x][y] = 1 if and only if x + y \\in \\Delta,

    where \\Delta is defined as

    \\Delta = \\{ (a, b), DDT_f[a][b] == 2  \\} ~.

    """
    n = int(log(len(f), 2))
    dim = 2**(2*n)
    d = ddt(f)
    gamma = [(a << n) | b
             for a, b in itertools.product(range(1, 2**n), range(0, 2**n))
             if d[a][b] == 2
    ]
    mat_content = []
    for x in range(0, dim):
        row = [0 for j in range(0, dim)]
        for y in gamma:
            row[oplus(x, y)] = 1
        mat_content.append(row)
    mat_gf2 = Matrix(GF(2), dim, dim, mat_content)
    return mat_gf2.rank()



# !SUBSECTION! Thickness related 

def thickness(basis, N):
    """Returns the thickness of the vector space with basis `basis`, where
    this vector space is a subspace of $\\F_2^{N+M}$ for the given `N` and
    for some $M$.

    The thickness is the rank of the projection of this space on its
    `N` bits of lowest weight.

    """
    MASK_N = sum(int(1 << i) for i in range(0, N))
    proj = [b & MASK_N for b in basis]
    return rank_of_vector_set(proj)


def thickness_spectrum(s, spaces=None):
    """Returns a dictionary containing the thickness spectra of the
    function whose LUT is the list `s`.

    If the spaces in the Walsh zeroes have already been extracted then
    it is possible to avoid their re-computation by passing them via
    the `spaces` input of this function.

    """
    N = int(log(len(s), 2))
    if spaces == None:
        spaces = get_lat_zeroes_spaces(s)
    result = defaultdict(int)
    for V in spaces:
        result[thickness(V, N)] += 1
    return dict(result)

def get_lat_zeroes_spaces(s, n_threads=DEFAULT_N_THREADS):
    """Returns a list containing the basis of each vector space of
    dimension n contained in the LAT zeroes of `s`.

    """
    return get_lat_zeroes_spaces_fast(s,
                                      int(log(len(s), 2)),
                                      int(n_threads))


# !SUBSECTION! Using the ortho-derivative

def ortho_derivative_label(f):
    """Returns a string representation of the differential and extended
    Walsh spectra of the ortho-derivative of the function given.

    Can only be applied to quadratic APN functions (this is not
    verified).

    """
    o = ortho_derivative(f)
    return "d{}; w{}".format(
        pretty_spectrum(differential_spectrum(o)),
        pretty_spectrum(walsh_spectrum(o), absolute=True),
    )
    

# !SUBSECTION! TU projection and decomposition

def tu_projection(s, t):
    """Returns a pair T, U of lists of lists such that

    s[(j << t) | i] == (U[i][j] << t) | T[j][i]

    i.e. such that T and U for the TU projection of s in the sense of
    https://eprint.iacr.org/2018/713

    """
    N = int(log(len(s), 2))
    mask = sum(int(1 << j) for j in range(0, N-t))
    T = [
        [-1 for i in range(0, 2**t)]
        for j in range(0, 2**(N-t))
    ]
    U = [
        [-1 for j in range(0, 2**(N-t))]
        for i in range(0, 2**t)
    ]
    for j in range(0, 2**(N-t)):
        for i in range(0, 2**t):
            y = s[(j << t) | i]
            T[j][i] = y & mask
            U[i][j] = y >> t
    return T, U


# !SECTION! Linear/Affine Equivalence 


# !SUBSECTION! XOR equivalence 

def xor_equivalence(f, g):
    """Returns a pair [k0, k1] of integers such that, for all x:

    f[x] = g[x + k0] + k1,

    where "+" denotes the bitwise XOR.
    """
    N = int(log(len(f), 2))
    for k0 in range(0, 2**N):
        k1 = oplus(f[0], g[k0])
        good = True
        for x in range(1, 2**N):
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
    """Returns a 64-char string obtained by hashing the base 16
    representation of each entry of the lookup table f with SHA-256.

    """
    hf = sha256()
    for x in f:
        hf.update(hex(x).encode('utf-8'))
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
    for c in range(0, len(f)):
        f_c = le_class_representative([oplus(f[x], c) for x in range(0, len(f))])
        d = hash_sbox(f_c)
        table_f[d] = f_c
        table_c[d] = c
    rs = []
    a = -1
    b = -1    
    for c in range(0, len(f)):
        g_c = le_class_representative([g[oplus(x, c)] for x in range(0, len(f))])
        d = hash_sbox(g_c)
        if d in table_c.keys():
            a = c
            b = table_c[d]
            rs = g_c
            break
    if a == -1:
        return []
    l_f = linear_equivalence([oplus(f[x], b) for x in range(0, len(f))],
                             rs)
    A_f, B_f = l_f[0], l_f[1]
    l_g = linear_equivalence([g[oplus(x, a)] for x in range(0, len(f))],
                             rs)
    A_g, B_g = l_g[0], l_g[1]
    A = A_g.inverse() * A_f
    B = B_f * B_g.inverse()
    a = apply_bin_mat(a, A.inverse())
    return [A, a, B, b]



# !SECTION! CCZ-equivalence 

# !SUBSECTION! CCZ-equivalence to a permutation

def ea_equivalent_permutation_mappings(f, spaces=None):
    """Returns the list of all the linear functions L such that f(x) +
    L(x) is a permutation, where L is in fact the matrix
    representation of this function.

    """
    N = int(log(len(f), 2))
    mask = sum((1 << i) for i in range(0, N))
    if spaces == None:
        spaces = get_lat_zeroes_spaces(f)
    result = []
    for V in spaces:
        if thickness(V, N) == N:
            L_lut = [-1 for x in range(0, 2**N)]
            full_space = linear_span(V)
            for x in full_space:
                L_lut[x & mask] = x >> N
            if -1 in L_lut:
                raise Exception("Problem in EA-equivalent mapping")
            else:
                result.append(
                    linear_function_lut_to_matrix(L_lut).transpose()
                )
    return result


def ccz_equivalent_permutations(f, 
                                number="all permutations", 
                                spaces=None,
                                minimize_ea_classes=False):
    """Returns a list of permutations that are CCZ-equivalent to
    `f`. 

    The behaviour of the function depends on the value of `number`:

    - if it is set to "all permutations" then at least one permutation
      per affine-equivalence class of permutations is returned; but

    - if it is set to "just one" then the output contains at most one
      permutation, namely the first one found.

    If the list of the vector space bases is known, it is possible to
    avoid its recomputation by passing it via the `spaces` argument.


    """
    N = int(log(len(f), 2))
    mask = sum(int(1 << i) for i in range(0, N))
    graph_f = [(x << N) | f[x] for x in range(0, 2**N)]
    if spaces == None:
        spaces = get_lat_zeroes_spaces(f)
    spaces_by_dimensions = defaultdict(list)
    for b in spaces:
        t1 = rank_of_vector_set([v >> N for v in b])
        t2 = rank_of_vector_set([v & mask for v in b])
        spaces_by_dimensions[(t1 << N) | t2].append(b)
    for dim_pairs in itertools.product(spaces_by_dimensions.keys(), 
                                       spaces_by_dimensions.keys()):
        t1, t2 = dim_pairs[0] >> N, dim_pairs[0] & mask
        u1, u2 = dim_pairs[1] >> N, dim_pairs[1] & mask
        if (t1 + u1) >= N and (t2 + u2) >= N:
            for b0 in spaces_by_dimensions[dim_pairs[0]]:
                for b1 in spaces_by_dimensions[dim_pairs[1]]:
                    if rank_of_vector_set(b0 + b1) == 2*N:
                        L = Matrix(GF(2), 2*N, 2*N, [tobin(x, 2*N) for x in b0 + b1])
                        g = apply_mapping_to_graph(f, L)
                        if number == "just one":
                            yield g
                            return 
                        else:
                            yield g
                            if minimize_ea_classes:
                                break
    



# !SUBSECTION! Exploring a CCZ-class


def apply_mapping_to_graph(f, L):
    n = int(log(len(f), 2))
    mask = sum(int(1 << i) for i in range(0, n))
    graph_f = [(x << n) | f[x] for x in range(0, 2**n)]
    L_map = FastLinearMapping(L)
    graph_g = [L_map(word) for word in graph_f]
    g = [-1 for x in range(0, 2**n)]
    for word in graph_g:
        x, y = word >> n, word & mask
        g[x] = y
    if -1 in g:
        raise Exception("The mapping is L not admissible for f")
    else:
        return g
    

def ccz_equivalent_function(f, V):
    """Assuming that V is a vector space of dimension n contained in the
    Walsh zeroes of f, applies a linear permutation L to the codebook of f
    which is such that L^T({(0, x), x \\in F_2^n}) = V.

    """
    n = int(log(len(f), 2))
    mask = sum(int(1 << i) for i in range(0, n))
    L_map = FastLinearMapping(get_generating_matrix(V, 2*n).transpose())
    graph_f = [(x << n) | f[x] for x in range(0, 2**n)]
    graph_g = [L_map(word) for word in graph_f]
    g = [-1 for x in range(0, 2**n)]
    for word in graph_g:
        x, y = word >> n, word & mask
        g[x] = y
    if -1 in g:
        raise Exception("V was not contained in the Walsh zeroes of f!")
    else:
        return g
    
    

def enumerate_ea_classes(f):
    """Returns a list containing at least one function from each of the
    EA-classes constituting the CCZ-class of `f`.

    Note that several functions in the same EA-class may be
    returned. Solving this problem is in fact an open *research*
    problem.

    """
    N = int(log(len(f), 2))
    mask = sum(int(1 << i) for i in range(0, N))
    graph_f = [(x << N) | f[x] for x in range(0, 2**N)]
    bases = get_lat_zeroes_spaces(f)
    result = []
    for b in bases:
        L_map = FastLinearMapping(get_generating_matrix(b, 2*N).transpose())
        graph_g = [L_map(word) for word in graph_f]
        g = [-1 for x in range(0, 2**N)]
        for word in graph_g:
            x, y = word >> N, word & mask
            g[x] = y
        if -1 in g:
            raise Exception("permutation ill defined!")
        else:
            result.append(g)
    return result


def ea_classes_in_the_ccz_class_of(f, include_start=False):
    """Returns an iterable that, when iterated over, will yield at least
    one function from each of the EA-classes constituting the
    CCZ-class of `f`.

    Note that several functions in the same EA-class may be
    returned. Solving this problem is in fact an open *research*
    problem.

    If `include_start` is set to False then the ea class of f is
    hopefully not returned. More precisely, the spaces with thickness
    0 are not considered.

    """
    N = int(log(len(f), 2))
    mask = sum(int(1 << i) for i in range(0, N))
    graph_f = [(x << N) | f[x] for x in range(0, 2**N)]
    z = lat_zeroes(f)
    for b in vector_spaces_bases_iterator(z, N, 2*N):
        if include_start or thickness(b, 2*N) > 0:
            L_map = FastLinearMapping(get_generating_matrix(b, 2*N).transpose())
            graph_g = [L_map(word) for word in graph_f]
            g = [-1 for x in range(0, 2**N)]
            for word in graph_g:
                x, y = word >> N, word & mask
                g[x] = y
            if -1 in g:
                raise Exception("CCZ map is ill defined!")
            else:
                yield g
    


# !SECTION! Tests

def print_result(n_valid, n_tested):
    verdict = "[success]"
    if (n_valid != n_tested):
        verdict = "[FAIL]"
    print("{} success rate: {}/{} = {:.03f}".format(
        verdict,
        n_valid,
        n_tested,
        float(n_valid)/n_tested))

# !SUBSECTION! Linear equivalence 

def check_linear_equivalence(f, g, A, B):
    """Returns True if and only if f = B o g o A."""
    for x in range(0, 2**N):
        y = apply_bin_mat(g[apply_bin_mat(x, A)], B)
        if y != f[x]:
            return False
    return True


def test_le_equivalence(N, verbose=False):
    from timeit import default_timer
    false_negatives = 0
    false_positives = 0
    n_tested = 500
    for i in range(0, n_tested):
        # checking if linearly equivalent permutations are identified
        g = random_permutation(N)
        A = rand_linear_permutation(N)
        B = rand_linear_permutation(N)
        f = [apply_bin_mat(g[apply_bin_mat(x, A)], B) for x in range(0, 2**N)]
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
                        print("[success]     LE {:0.4f}s (other matrices found)".format(elapsed))
                else:
                    false_negatives += 1
                    if verbose:
                        print("[FAIL]     LE {:0.4f}s (wrong matrices found)".format(elapsed))
            else:
                if verbose:
                    print("[success]     LE {:0.4f}s".format(elapsed))
        else:
            false_negatives += 1
            if verbose:
                print("[FAIL]     LE {:0.4f} (nothing found)".format(
                    computation_end - computation_start))
        # checking if non-linearly equivalent functions are identified
        g = random_permutation(N)
        result = linear_equivalence(f, g)
        if len(result) == 0:
            if verbose:
                print("[success] non-LE {:0.4f}".format(elapsed))
        else:
            if check_linear_equivalence(f, g, result[0], result[1]):
                if verbose:
                    print("[success] act.LE {:0.4f}".format(elapsed))
            else:
                false_positives += 1
                if verbose:
                    print("[FAIL] matrices found for non-LE permutations")
    print("* testing if LE functions are identified correctly (with correct linear permutations)")
    print_result(n_tested-false_negatives, n_tested)
    print("* testing if NON-LE functions are identified correctly")
    print_result(n_tested-false_positives, n_tested)



# !SUBSECTION! LE representative 

def test_le_repr(N, verbose=False):
    import diff_lin
    from timeit import default_timer
    n_valid = 0
    n_tested = 50
    print("* Testing whether f and le_class_representative(f) are LE")
    for i in range(0, n_tested):
        f = random_permutation(N)
        computation_start = default_timer()
        g = le_class_representative(f)
        computation_end = default_timer()
        if len(linear_equivalence(f, g)) != 0:
            n_valid += 1
            if verbose:
                print("[success] {:0.4f}s".format(computation_end - computation_start))
        else:
            print("[FAIL]")
    print_result(n_valid, n_tested)
    print("* testing whether two linear equivalent functions have the same representative")
    n_valid = 0
    for i in range(0, n_tested):
        f = random_permutation(N)
        A = rand_linear_permutation(N)
        B = rand_linear_permutation(N)
        g = [apply_bin_mat(f[apply_bin_mat(x, A)], B) for x in range(0, 2**N)]
        rs_f = le_class_representative(f)
        rs_g = le_class_representative(g)
        identical = True
        for x in range(0, 2**N):
            if rs_f[x] != rs_g[x]:
                identical = False
                break
        if identical:
            n_valid += 1
            if verbose:
                print("[success]")
        else:
            if verbose:
                print("[FAIL] representatives don't match")
                print(rs_f, pretty_spectrum(diff_lin.differential_spectrum(rs_f)))
                print(rs_g, pretty_spectrum(diff_lin.differential_spectrum(rs_g)))
    print_result(n_valid, n_tested)
    

# !SUBSECTION! Affine equivalence
    
def check_affine_equivalence(f, g, A, a, B, b):
    """Checks whether f(x) = (B o g o A)(x + a) + b"""
    for x in range(0, 2**N):
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
    for i in range(0, n_tested):
        # checking if linearly equivalent permutations are identified
        f = random_permutation(N)
        A = rand_linear_permutation(N)
        a = randint(0, 2**N-1)
        B = rand_linear_permutation(N)
        b = randint(0, 2**N-1)
        g = [oplus(apply_bin_mat(f[apply_bin_mat(oplus(x, a), A)], B), b)
             for x in range(0, 2**N)]
        computation_start = default_timer()
        result = affine_equivalence(f, g)
        computation_end = default_timer()
        elapsed = computation_end - computation_start
        if len(result) > 1:
            if not check_affine_equivalence(f, g, result[0], result[1], result[2], result[3]):
                false_negatives += 1
                if verbose:
                    print("[FAIL] wrong affine permutations")
            else:
                if verbose:
                    print("[success]     AE {:0.4f}".format(elapsed))

        else:
            false_negatives += 1
            if verbose:
                print("[FAIL]     AE {:0.4f}s (nothing found)".format(elapsed))
        # checking if non-affine equivalent functions are identified
        g = random_permutation(N)
        result = affine_equivalence(f, g)
        if len(result) == 0:
            if verbose:
                print("[success] non-AE {:0.4f}s".format(elapsed))
        else:
            if check_affine_equivalence(f, g, result[0], result[1], result[2], result[3]):
                if verbose:
                    print("[success] act.AE {:0.4f}".format(elapsed))
            else:
                false_positives += 1
                if verbose:
                    print("[FAIL] matrices found for non-LE permutations")
    print("* testing if AE functions are identified correctly (with correct affine permutations)")
    print_result(n_tested-false_negatives, n_tested)
    print("* testing if NON-LE functions are identified correctly")
    print_result(n_tested-false_positives, n_tested)



# !SUBSECTION!  Test CCZ-equivalent permutation

def test_ea_permutations():
    for N in [4, 5]:
        F = GF(2**N, name="a")
        inv = [(F.fetch_int(x)**(2**N-2)).integer_representation()
               for x in range(0, 2**N)]
        print("== " + str(N))
        for L in ea_equivalent_permutation_mappings(inv):
            print(L.str() + "\n")

def test_ccz_permutations(number="all permutations"):
    N = 6
    F = GF(2**N, name="a")
    # generating the Kim mapping
    kim = []
    for x_i in range(0, 2**N):
        x = F.fetch_int(x_i)
        y = x**3 + x**10 + F.gen()*x**24
        kim.append(y.integer_representation())
    permutations = ccz_equivalent_permutations(kim, number=number)
    for i, p in enumerate(permutations):
        print("{:2d} {} {} {}".format(
            i,
            is_permutation(p),
            pretty_spectrum(differential_spectrum(p)),
            pretty_vector(p)            
        ))
    print("total: {}".format(len(permutations)))

def test_enumerate_ea():
    N = 8
    F = GF(2**N, name="a")
    # generating the Kim mapping
    kim = []
    for x_i in range(0, 2**N):
        x = F.fetch_int(x_i)
        y = x**3 + x**10 + F.gen()*x**24
        kim.append(y.integer_representation())
    classes = enumerate_ea_classes(kim)
    for f in classes:
        print(str(algebraic_degree(f)) + pretty_spectrum(thickness_spectrum(f)))
    print("total: " + str(len(classes)))
    

def test_ea_classes():
    N = 8
    F = GF(2**N, name="a")
    # generating the Kim mapping
    kim = []
    for x_i in range(0, 2**N):
        x = F.fetch_int(x_i)
        y = x**3 + x**10 + F.gen()*x**24
        kim.append(y.integer_representation())

    total = 0
    for f in ea_classes_in_the_ccz_class_of(kim):
        print(str(algebraic_degree(f)) + pretty_spectrum(thickness_spectrum(f)))
        total += 1
    print("total: " + str(total))
    

# !SECTION! Running tests

if __name__ == '__main__':
    # test_ea_permutations()
    # test_ccz_permutations(number="just one")
    # test_enumerate_ea()
    test_ea_classes()
    
    # import sys
    # N = int(sys.argv[1])
    # print("=== Linear Equivalence ===")
    # test_le_equivalence(N, verbose=True)
    # print("\n=== Linear Representative ===")
    # test_le_repr(N, verbose=False)
    # print("\n=== Affine Equivalence ===")
    # test_ae_equivalence(N, verbose=True)
