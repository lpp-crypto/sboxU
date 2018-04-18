#!/usr/bin/sage

from sage.all import Matrix, GF, vector, log, randint, Graph, ceil
from sage.graphs.cliquer import all_max_clique

import itertools
import random
from sboxu_cpp import *
from display import *
from linear import *

DEFAULT_N_THREADS  = 16


def inverse(s):
    result = [0 for i in xrange(0, len(s))]
    for x in xrange(0, len(s)):
        result[s[x]] = x
    return result
    

# !SECTION! Getting zeroes

def proj_lat_zeroes(s):
    result = []
    for b in xrange(1, len(s)):
        w = fourier_transform([scal_prod(b, s[x]) for x in xrange(0, len(s))])
        for c in w:
            if c == 0:
                result.append(b)
                break
    return result
                
# !SECTION! Finding vector spaces of zeroes


def indicator_function(l):
    result = defaultdict(int)
    for x in l:
        result[x] = 1
    return result


def intersection(s1, s2):
    d2 = defaultdict(int)
    for x in s2:
        d2[x] = 1
    result = []
    for x in s1:
        if d2[x] == 1:
            result.append(x)
    return result


def list_to_integer(l, N):
    l.sort()
    result = 0
    for x in l:
        result = (result << N) | x
    return result
    

def integer_to_list(x, N):
    result = []
    y = x
    mask = sum(1 << i for i in xrange(0, N))
    while y != 0:
        result.append(y & mask)
        y = y >> N
    result.sort()
    return result
    


def extract_bases(z, dimension, n_threads=DEFAULT_N_THREADS):
    return extract_bases_fast(z, int(dimension), n_threads)


def extract_ae_bases(l):
    result = []
    N = int(log(len(l), 2))
    for a in xrange(0, len(l)):
        if l[a][1] == 0:
            result += extract_ae_bases_rec(l, N, [a], [[0, 0], [a, 1]])
    return result


def extract_ae_bases_rec(l, N, base, span_base):
    """l is the LAT of a function, base is the basis currently considered and
    span_base is its span."""
    if len(base) == N:
        return [base]
    else:
        b = (1 << len(base))
        result = []
        for a in xrange(0, len(l)):
            if l[a][b] == 0:
                new_span_base = copy(span_base)
                valid_guess = True
                for v in span_base:
                    a_prime = oplus(v[0], a)
                    b_prime = oplus(v[1], b)
                    if l[a_prime][b_prime] != 0:
                        valid_guess = False
                        break
                    new_span_base.append([a_prime, b_prime])
                if valid_guess:
                    # print new_span_base
                    result += extract_ae_bases_rec(l, N, base + [a], new_span_base)
        return result


def minimizing_offset(b, N):
    s = span(b)
    best_offset = 0
    min_offset_s = s
    for offset in xrange(0, 2**N):
        new_s = [oplus(x, offset) for x in s]
        if new_s < min_offset_s:
            min_offset_s = new_s
            best_offset = offset
    return best_offset
        


def offsets_rec(N, z_basis, offsets, z_columns, z_columns_charac, current_intersection, n_threads=DEFAULT_N_THREADS):
    d = len(z_basis)
    if len(current_intersection) < 2**(N-d):
        return []
    if len(offsets) == N-d:
        return [[offsets, current_intersection]]
    index = len(offsets)
    b = z_basis[index]

    result = []
    for a in z_columns[b]:
        new_intersection_indicator = indicator_function(current_intersection)
        for i in xrange(0, 2**index):
            new_offset = a
            new_b = b
            for j in xrange(0, index):
                if (i >> j) & 1 == 1:
                    new_b = oplus(new_b, z_basis[j])
                    new_offset = oplus(new_offset, offsets[j])
            for y in current_intersection:
                if z_columns_charac[new_b][oplus(y, new_offset)] != 1:
                    new_intersection_indicator[y] = 0
        new_intersection = [y for y in new_intersection_indicator.keys()
                            if new_intersection_indicator[y] == 1]
        if len(new_intersection) >= 2**(N-d)-1:
            result += offsets_rec(N,
                                  z_basis,
                                  offsets + [a],
                                  z_columns,
                                  z_columns_charac,
                                  new_intersection,
                                  n_threads=n_threads)
            
    return result
                

def find_permutation_spaces_quadratic(l, N, n_threads=DEFAULT_N_THREADS):
    z_0_indicator = defaultdict(int)
    for a in xrange(0, 2**N):
        for b in xrange(0, 2**N):
            if l[a][b] == 0:
                z_0_indicator[b] = 1
    z_0 = [b for b in xrange(0, 2**N) if z_0_indicator[b] == 1]
    spaces = extract_bases(z_0, int(N/2), n_threads=n_threads)
    return spaces
    

def find_permutation_spaces(l, N, n_threads=16):
    z_0 = [b for b in xrange(0, 2**N) if l[0][b] % 2**N == 0]
    z_0_bases = []
    for d in xrange(int(N/2), N-1):
        new_bases = extract_bases(z_0, d, n_threads=n_threads)
        if len(new_bases) == 0:
            break
        else:
            z_0_bases += new_bases

    big_enough_space_found = False
    for base in z_0_bases:
        if len(base) >= int(N/2):
            big_enough_space_found = True
    if not big_enough_space_found:
        return []

    z_columns = [[] for b in xrange(0, 2**N)]
    for b in xrange(0, 2**N):
        for a in xrange(0, 2**N):
            if l[a][b] % 2**N == 0:
                v = (a << N) | b
                z_columns[b].append(a)
    z_columns_charac = [indicator_function(z_columns[b]) for b in xrange(0, 2**N)]
                
    for b1, b2 in itertools.combinations(z_0_bases, 2):
        if len(b1) + len(b2) == N and rank_of_vector_set(b1 + b2, 2*N) == N:
            print pretty_vector(b1), pretty_vector(b2)
            
            correction_result1 = offsets_rec(N,
                                             b1,
                                             [],
                                             z_columns,
                                             z_columns_charac,
                                             range(0, 2**N),
                                             n_threads=n_threads)
            print len(correction_result1)
            correction_result2 = offsets_rec(N,
                                             b2,
                                             [],
                                             z_columns,
                                             z_columns_charac,
                                             range(0, 2**N),
                                             n_threads=n_threads)
            print len(correction_result2)
            # if len(correction_result) > 0:
            #     correction, inter = correction_result
            #     print correction
            #     for base_in_intersection in extract_bases(inter, N-d, n_threads=n_threads):
            #         space = []
            #         for b_index in xrange(0, 2**d):
            #             b = 0
            #             cor = 0
            #             for j in xrange(0, d):
            #                 if (b_index >> j) & 1 == 1:
            #                     b = oplus(b, proj_base[j])
            #                     cor = oplus(cor, correction[j])
            #             for a_index in xrange(0, 2**(N-d)):
            #                 a = cor
            #                 for j in xrange(0, N-d):
            #                     if (a_index >> j) & 1 == 1:
            #                         a = oplus(a, base_in_intersection[j])
            #                 space.append((a << N) | b)
            #         bases = extract_bases(space, N, n_threads=16)
            #         candidate_bases += bases

    # print len(candidate_bases)
    # mask = sum(1 << i for i in xrange(0, N))
    # for c in candidate_bases:
    #     counter = 0
    #     for x in span(c):
    #         a = x >> N
    #         b = x & mask
    #         if l[a][b] % 2**N == 0:
    #             counter += 1
    #     print [(hex(x >> N), hex(x & mask)) for x in c], rank_of_vector_set(c, 2*N), counter
    # result = []
    # for v1, v2 in itertools.product(candidate_bases, candidate_bases):
    #     if rank_of_vector_set(v1 + v2, 2*N) == 2*N:
    #         result.append([v1, v2])
    #         # print pretty_vector(v1, template="{:x}"), pretty_vector(v2, template="{:x}")
    # return result
    return []


def extract_direct_sum(list_of_spaces, N):
    space_pairs = []
    for va, vb in itertools.combinations(list_of_spaces, 2):
        # basic tests
        good = True
        if len(va[0])*len(vb[0]) < 2**N:
            good = False
        if good:
            if len(va[1])*len(vb[1]) < 2**N:
                good = False
        if good:
            if rank_of_vector_set(va[0] + vb[0], N) != N:
                good = False
        if good:
            if rank_of_vector_set(va[1] + vb[1], N) != N:
                good = False
        if good:
            # basic tests passed
            # more complicated cases
            # -- the y coordinate is good but the x coordinate is too big
            if len(va[1])*len(vb[1]) == 2**N and len(va[0])*len(vb[0]) > 2**N:
                if len(vb[0]) > len(va[0]): # case where vb needs shringking
                    # print "1", va[0], vb[0]
                    new_vb0_indicator = defaultdict(int)
                    for x in vb[0]:
                        new_vb0_indicator[x] = 1
                    for v in va[0]:
                        if v != 0:
                            for x in vb[0]:
                                new_vb0_indicator[oplus(x, v)] = 0
                    new_b0 = [x for x in new_vb0_indicator.keys() if new_vb0_indicator[x] == 1]
                    vb[0] = new_b0
                elif len(va[0]) > len(vb[0]): # case where va needs shringking
                    # print "2"
                    new_va0_indicator = defaultdict(int)
                    for x in va[0]:
                        new_va0_indicator[x] = 1
                    for v in vb[0]:
                        if v != 0:
                            for x in va[0]:
                                new_va0_indicator[oplus(x, v)] = 0
                    new_a0 = [x for x in new_va0_indicator.keys() if new_va0_indicator[x] == 1]
                    va[0] = new_a0
            # -- the x coordinate is good but the y coordinate is too big
            if len(va[0])*len(vb[0]) == 2**N and len(va[1])*len(vb[1]) > 2**N:
                if len(vb[1]) > len(va[1]): # case where vb needs shringking
                    # print "3"
                    new_vb1_indicator = defaultdict(int)
                    for x in vb[1]:
                        new_vb1_indicator[x] = 1
                    for v in va[1]:
                        if v != 0:
                            for x in vb[1]:
                                new_vb1_indicator[oplus(x, v)] = 0
                    new_b1 = [x for x in new_vb1_indicator.keys() if new_vb1_indicator[x] == 1]
                    vb[1] = new_b1
                elif len(va[1]) > len(vb[1]): # case where va needs shringking
                    # print "4"
                    new_va1_indicator = defaultdict(int)
                    for x in va[1]:
                        new_va1_indicator[x] = 1
                    for v in vb[1]:
                        if v != 0:
                            for x in va[1]:
                                new_va1_indicator[oplus(x, v)] = 0
                    new_a1 = [x for x in new_va1_indicator.keys() if new_va1_indicator[x] == 1]
                    va[1] = new_a1 
            # simple case: all spaces are just the right size
            if len(va[0])*len(vb[0]) == 2**N and len(va[1])*len(vb[1]) == 2**N:
                space_pairs.append([va, vb])
    return space_pairs


# !SECTION! Exploiting vector spaces 


# !SUBSECTION! Obtaining T and U (when no linear layers are involved) 

def get_tu_open(s, t):
    """Assumes that s is the codebook of a lopsided open butterfly with
    branch widths of t and N-t, and returns the corresponding mini-block
    ciphers T and U.

    """
    n = int(log(len(s), 2))
    mask_t = sum(int(1 << i) for i in xrange(0, t))
    T = [[0 for l in xrange(0, 2**t)] for r in xrange(0, 2**(n-t))]
    U = [[0 for r in xrange(0, 2**(n-t))] for l in xrange(0, 2**t)]
    for l in xrange(0, 2**t):
        for r in xrange(0, 2**(n-t)):
            x = (l << (n-t)) | r
            # x = (r << t) | l
            y_l, y_r = s[x] >> t, s[x] & mask_t
            T[r][l] = y_r
            U[y_r][r] = y_l
    return T, U


def get_tu_closed(s, t):
    """Assumes that s is the codebook of a lopsided closed butterfly with
    branch widths of t and N-t, and returns the corresponding
    mini-block ciphers T and U.

    """
    n = int(log(len(s), 2))
    mask_t = sum(int(1 << i) for i in xrange(0, t))
    T = [[0 for l in xrange(0, 2**t)] for r in xrange(0, 2**(n-t))]
    U = [[0 for r in xrange(0, 2**(n-t))] for l in xrange(0, 2**t)]
    for l in xrange(0, 2**t):
        for r in xrange(0, 2**(n-t)):
            x = (l << (n-t)) | r
            y_l, y_r = s[x] >> t, s[x] & mask_t
            T[r][l] = y_r
            U[l][r] = y_l
    return T, U


# !SUBSECTION! Exploiting spaces which are cartesian products

def tu_decomposition(s, v, verbose=False):
    """Using the knowledge that v is a subspace of Z_s of dimension n, a
    TU-decomposition (as defined in [Perrin17]) of s is performed and T
    and U are returned.

    """
    N = int(log(len(s), 2))
    # building B
    basis = extract_basis(v[0], N)
    t = len(basis)
    basis = complete_basis(basis, N)
    B = Matrix(GF(2), N, N, [tobin(x, N) for x in reversed(basis)])
    if verbose:
        print "B=  (rank={})\n{}".format(B.rank(), B.str())
    # building A
    basis = complete_basis(extract_basis(v[1], N) , N)
    A = Matrix(GF(2), N, N, [tobin(x, N) for x in reversed(basis)])
    if verbose:
        print "A=  (rank={})\n{}".format(A.rank(), A.str())
    # building linear equivalent s_prime
    s_prime = [apply_bin_mat(s[apply_bin_mat(x, A.inverse())], B)
               for x in xrange(0, 2**N)]
    # TU decomposition of s'
    T, U = get_tu_open(s_prime, t)
    if verbose:
        print "T=["
        for i in xrange(0, 2**(N-t)):
            print "  {} {}".format(T[i], is_permutation(T[i]))
        print "]\nU=["
        for i in xrange(0, 2**t):
            print "  {} {}".format(U[i], is_permutation(U[i]))
        print "]"
    return T, U


def get_ccz_equivalent_permutation_cartesian(s, v0, v1, verbose=False):
    """Takes as input two vector spaces v0 and v1 of the set of zeroes in
    the LAT of s, each being the cartesian product of two spaces, and
    returns a permutation CCZ equivalent to s obtained using these CCZ
    spaces.

    """
    N = int(log(len(s), 2))
    # building B
    e1 = extract_basis(v0[0], N)
    t = len(e1)
    e2 = extract_basis(v1[0], N)
    B = Matrix(GF(2), N, N, [tobin(x, N) for x in e1+e2])
    if verbose:
        print "B=  (rank={})\n{}".format(B.rank(), B.str())
    # building A
    e1 = extract_basis(v0[1], N) 
    e2 = extract_basis(v1[1], N)
    A = Matrix(GF(2), N, N, [tobin(x, N) for x in e1+e2])
    if verbose:
        print "A=  (rank={})\n{}".format(A.rank(), A.str())
    # building linear equivalent s_prime
    s_prime = [apply_bin_mat(s[apply_bin_mat(x, A.inverse())], B)
               for x in xrange(0, 2**N)]
    # TU decomposition of s'
    T_inv, U = get_tu_closed(s_prime, t)
    if verbose:
        print "T_inv=["
        for i in xrange(0, 2**t):
            print "  {} {}".format(T_inv[i], is_permutation(T_inv[i]))
        print "]\nU=["
        for i in xrange(0, 2**t):
            print "  {} {}".format(U[i], is_permutation(U[i]))
        print "]"
    # TU-unfolding
    T = [inverse(row) for row in T_inv]
    result = [0 for x in xrange(0, 2**N)]
    for l in xrange(0, 2**t):
        for r in xrange(0, 2**(N-t)):
            x = (l << (N-t)) | r
            y_r = T[r][l]
            y_l = U[y_r][r]
            result[x] = (y_l << (t)) | y_r
    return result


def get_ccz_equivalent_function_cartesian(s, v, verbose=False):
    """Takes as input a vector space v of the set of zeroes in
    the LAT of s, v being the cartesian product of two spaces, and
    returns a function CCZ equivalent to s obtained using this CCZ
    space.

    """
    N = int(log(len(s), 2))
    T_inv, U = tu_decomposition(s, v, verbose=verbose)
    t = int(log(len(T_inv[0]), 2))
    T = [inverse(row) for row in T_inv]
    result = [0 for x in xrange(0, 2**N)]
    for l in xrange(0, 2**t):
        for r in xrange(0, 2**(N-t)):
            x = (l << (N-t)) | r
            y_r = T[r][l]
            y_l = U[l][r]
            result[x] = (y_l << t) | y_r
    return result


# !SUBSECTION! Exploiting all spaces

def get_ccz_equivalent_function(l, v, verbose=False):
    N = int(log(len(l), 2))
    mask = sum(int(1 << i) for i in xrange(0, N))
    L = complete_basis(v, 2*N)
    L = Matrix(GF(2), 2*N, 2*N, [tobin(x, 2*N) for x in L]).transpose()
    if verbose:
        print "\n\nrank={}\n{}".format(L.rank(), L.str())
    new_l = [[0 for b in xrange(0, 2**N)] for a in xrange(0, 2**N)]
    for a, b in itertools.product(xrange(0, 2**N), xrange(0, 2**N)):
        v = apply_bin_mat((a << N) | b, L)
        a_prime, b_prime = v >> N, v & mask
        new_l[a][b] = l[a_prime][b_prime]
    return invert_lat(new_l)



def get_ccz_equivalent_permutation(l, v0, v1, verbose=False):
    N = int(log(len(l), 2))
    mask = sum(int(1 << i) for i in xrange(0, N))
    L = v0 + v1
    L = Matrix(GF(2), 2*N, 2*N, [tobin(x, 2*N) for x in L]).transpose()
    if verbose:
        print "\n\nrank={}\n{}".format(L.rank(), L.str())
    new_l = [[0 for b in xrange(0, 2**N)] for a in xrange(0, 2**N)]
    for a, b in itertools.product(xrange(0, 2**N), xrange(0, 2**N)):
        v = apply_bin_mat((a << N) | b, L)
        a_prime, b_prime = v >> N, v & mask
        new_l[a][b] = l[a_prime][b_prime]
    return invert_lat(new_l)


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
