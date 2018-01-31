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
    

# !SECTION! Finding vector spaces of zeroes

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
    

def lat_zeroes(s):
    l = lat(s)
    words = [0]
    for a, b in itertools.product(xrange(0, 2**N), xrange(0, 2**N)):
        if l[a][b] == 0:
            words.append((a << N) | b)
    return words

def extract_bases(z, dimension, n_threads=DEFAULT_N_THREADS):
    return extract_bases_fast(z, dimension, n_threads)

# def clique_based_filter(z):
#     g = Graph()
#     g.add_vertices(z)
#     indicator_z = defaultdict(int)
#     for x in z:
#         indicator_z[x] = 1
#     for a in z:
#         g.add_edge(a, a)
#     for a, b in itertools.combinations(z, r=2):
#         u = oplus(a, b)
#         if indicator_z[u] == 1:
#             g.add_edge(a, b)
#     c = all_max_clique(g)
#     return c


# def indicator_function(z):
#     result = defaultdict(int)
#     for x in z:
#         result[x] = 1
#     return result


# def extract_vector(z, a):
#     if len(z) == 0:
#         return []
#     indicator_z = indicator_function(z)
#     result = []
#     for x in z:
#         y = oplus(x, a)
#         if y > x and indicator_z[y] == 1:
#             result.append(x)
#     return result


# def extract_spaces(z, d):
#     if d < 1:
#         return []
#     result = []
#     for a in z:
#         z_a = extract_vector(z, a)
#         filtered_z_a = [x for x in z_a if x > a]
#         result += extract_spaces_rec([a], filtered_z_a, d)
#     return result


# def extract_spaces_rec(basis, z, d):
#     result = []
#     if len(basis) == 0:
#         raise Exception("extract_spaces_rec does not allow an empty basis.")
#     if len(z) == 0 or len(z) < 2**(d - len(basis)) - 1:
#         return []
#     if len(basis) == d-1:
#         for a in z:
#             if a > basis[-1]:
#                 result.append(basis + [a])
#     else:
#         for a in z:
#             if a > basis[-1]:
#                 z_a = extract_vector(z, a)
#                 filtered_z_a = [x for x in z_a if x > a]
#                 result += extract_spaces_rec(basis + [a], filtered_z_a, d)
#     return result
        
    
def extract_zero_spaces(s, verbose=False, transpose=True):
    """Computes the LAT of s and returns a list of vector spaces S_i where
    S_i is a subspaces of $[0, 2^N] \times [0, 2^N]$ of dimension n.

    """
    N = int(log(len(s), 2))
    mask_N = sum(int(1 << i) for i in xrange(0, N))
    
    l = lat(s)
    
    search_space = ["lat"]
    if transpose:
        search_space.append("lat transpose")

    tu_spaces_all = []
    for search in search_space:
        if verbose:
            print "\n* search in " + search
            print "** Looking for vector spaces of zeroes in rows..."
            
        row_zero_spaces = defaultdict(list)
        row_zero_cliques = []
        for a in xrange(0, 2**N):
            row_zeroes = []
            for b in xrange(0, 2**N):
                if search != "lat":
                    if l[b][a] % 2**N == 0:
                        row_zeroes.append(b)
                else:
                    if l[a][b] % 2**N == 0:
                        row_zeroes.append(b)
            cliques = clique_based_filter(row_zeroes)
            big_cliques = []
            for c in cliques:
                if len(c) >= 2**(int(N/2)):# and 0 in c:
                    big_cliques.append(c)
            row_zero_cliques.append(big_cliques)
            if verbose:
                print a, big_cliques

        if verbose:
            print "  ... [DONE]"
            print "** Pairing rows to find common (large enough) vector spaces of zeroes..."
            
        row_zero_spaces = defaultdict(list)
        for i, j in itertools.combinations(xrange(0, 2**N), 2):
            k = oplus(i, j)
            if k > i and k > j:
                for c1, c2 in itertools.product(row_zero_cliques[i], row_zero_cliques[j]):
                    inter_ij = intersection(c1, c2)
                    if len(inter_ij) >= 2**(int(N/2)):
                        for c3 in row_zero_cliques[k]:
                            inter_ijk = intersection(inter_ij, c3)
                            if len(inter_ijk) >= 2**(int(N/2)):
                                inter_ijk.sort()
                                space_identifier = sum(inter_ijk[i_b] << (i_b * N) for i_b in xrange(0, len(inter_ijk)))
                                if verbose:
                                    print "  {:02x} {:02x} {:02x}   {}".format(
                                        i, j, k,
                                        pretty_vector(inter_ijk),
                                    )
                                for v in [0, i, j, k]:
                                    if v not in row_zero_spaces[space_identifier]:
                                        row_zero_spaces[space_identifier].append(v)
        
        if verbose:
            print "  ... [DONE]"
            print "** Looking for (large enough) vector spaces in row indices sharing same spaces of zeroes ..."

        tu_spaces = []
        for space_identifier in row_zero_spaces.keys():
            original_space = []
            y = space_identifier
            while y != 0:
                original_space.append(int(y & mask_N))
                y = y >> N
            original_space.sort()
            all_cliques = clique_based_filter(row_zero_spaces[space_identifier])
            for c in all_cliques:
                if len(c) * len(original_space) >= 2**N:
                    if search == "lat":
                        tu_spaces.append([original_space, c])
                    else:
                        tu_spaces.append([c, original_space])
                        
        for ss in tu_spaces:
            if ss not in tu_spaces_all:
                tu_spaces_all.append(ss)
                if verbose:
                    print "  {}, {}".format(
                        pretty_vector(ss[0], width=ceil(N/4)),
                        pretty_vector(ss[1], width=ceil(N/4)),
                    )
                
        if verbose:
            print "... [DONE]"

    return tu_spaces_all


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
