from sage.all import *
from sage.crypto import sboxes
from sboxU import *
import random
import itertools
import time

banff_list = [[0, 1, 8, 15, 27, 14, 35, 48, 53, 39, 43, 63, 47, 41, 1, 1, 41, 15, 15, 47, 52, 6, 34, 22, 20, 33, 36, 23, 8, 41, 8, 47, 36, 52, 35, 53, 35, 39, 20, 22, 33, 34, 48, 53, 39, 48, 6, 23, 22, 33, 63, 14, 23, 52, 14, 43, 27, 63, 36, 6, 27, 43, 20, 34], [0, 58, 31, 63, 13, 7, 44, 60, 59, 17, 2, 50, 61, 39, 58, 58, 39, 63, 63, 61, 30, 54, 56, 10, 45, 37, 19, 1, 31, 39, 31, 61, 19, 30, 44, 59, 44, 17, 45, 10, 37, 56, 60, 59, 17, 60, 54, 1, 10, 37, 50, 7, 1, 30, 7, 2, 13, 50, 19, 54, 13, 2, 45, 56], [0, 17, 60, 37, 20, 46, 0, 50, 6, 53, 8, 51, 40, 48, 14, 30, 0, 10, 48, 50, 2, 35, 26, 51, 61, 21, 63, 31, 5, 6, 47, 36, 48, 27, 3, 32, 5, 5, 30, 22, 15, 6, 14, 15, 0, 34, 41, 3, 26, 42, 37, 29, 57, 34, 46, 61, 30, 12, 19, 9, 7, 62, 34, 19], [0, 53, 39, 35, 55, 58, 8, 52, 23, 14, 56, 16, 31, 62, 40, 56, 54, 7, 50, 50, 50, 59, 46, 22, 59, 38, 55, 27, 0, 37, 20, 0, 11, 60, 26, 28, 56, 55, 49, 15, 58, 33, 35, 9, 54, 21, 55, 37, 33, 18, 19, 17, 33, 42, 11, 49, 10, 21, 48, 30, 53, 18, 23, 1], [0, 2, 52, 16, 25, 21, 54, 28, 37, 12, 35, 44, 16, 55, 13, 12, 25, 2, 1, 60, 28, 9, 31, 44, 53, 5, 31, 9, 28, 34, 45, 53, 34, 40, 29, 49, 21, 17, 49, 19, 26, 59, 23, 16, 1, 46, 23, 30, 5, 22, 22, 35, 46, 51, 38, 29, 52, 12, 21, 11, 51, 5, 9, 25], [0, 1, 9, 52, 6, 1, 13, 54, 42, 16, 39, 33, 14, 50, 1, 1, 0, 31, 60, 31, 58, 35, 4, 33, 49, 21, 9, 17, 41, 11, 19, 13, 16, 49, 14, 19, 28, 59, 0, 27, 29, 7, 7, 33, 51, 47, 43, 11, 1, 62, 42, 41, 49, 8, 24, 29, 23, 19, 56, 0, 5, 7, 40, 22], [0, 8, 50, 42, 4, 0, 24, 12, 3, 39, 56, 12, 33, 9, 52, 12, 59, 56, 41, 58, 5, 10, 57, 38, 11, 36, 16, 47, 19, 48, 38, 21, 8, 50, 40, 2, 63, 9, 49, 23, 35, 53, 10, 12, 50, 40, 53, 63, 63, 14, 63, 30, 50, 15, 28, 49, 39, 58, 46, 35, 12, 29, 43, 42], [0, 50, 59, 44, 53, 14, 33, 63, 61, 28, 31, 27, 34, 10, 47, 34, 2, 30, 21, 44, 61, 40, 5, 53, 58, 53, 52, 30, 47, 41, 14, 45, 3, 11, 62, 19, 50, 51, 32, 4, 28, 7, 56, 6, 7, 21, 12, 59, 44, 10, 61, 62, 23, 56, 41, 35, 54, 3, 62, 46, 39, 27, 0, 25], [0, 1, 41, 31, 37, 52, 22, 48, 28, 39, 29, 17, 32, 11, 59, 39, 19, 27, 41, 22, 43, 51, 11, 36, 38, 20, 52, 49, 7, 37, 15, 26, 52, 39, 10, 46, 26, 25, 62, 10, 4, 45, 18, 12, 51, 10, 63, 49, 22, 12, 59, 22, 37, 47, 18, 47, 15, 47, 10, 29, 37, 21, 58, 61], [0, 57, 9, 18, 22, 34, 39, 49, 29, 28, 30, 61, 34, 46, 25, 55, 5, 44, 47, 36, 29, 57, 15, 9, 28, 13, 60, 15, 45, 49, 53, 11, 7, 9, 0, 44, 1, 2, 62, 31, 62, 8, 51, 39, 17, 42, 36, 61, 3, 29, 39, 27, 11, 24, 23, 38, 62, 24, 16, 20, 31, 52, 9, 0], [0, 52, 13, 5, 9, 20, 52, 21, 45, 40, 38, 31, 39, 11, 28, 12, 21, 9, 9, 41, 20, 33, 56, 49, 53, 24, 47, 62, 55, 51, 29, 37, 45, 26, 14, 5, 18, 12, 1, 35, 44, 42, 9, 51, 16, 63, 5, 22, 27, 4, 41, 10, 44, 26, 46, 36, 23, 57, 35, 49, 35, 36, 39, 28], [0, 45, 13, 10, 51, 58, 39, 4, 10, 48, 33, 49, 62, 32, 12, 56, 29, 37, 7, 21, 62, 34, 61, 11, 8, 39, 52, 49, 44, 39, 9, 40, 18, 37, 36, 57, 45, 62, 2, 59, 63, 31, 47, 37, 7, 3, 14, 32, 30, 60, 63, 55, 49, 55, 9, 37, 44, 25, 43, 52, 4, 21, 26, 33], [0, 17, 28, 34, 59, 2, 20, 2, 1, 20, 8, 50, 56, 5, 2, 16, 35, 13, 50, 51, 11, 13, 41, 0, 50, 24, 54, 51, 24, 26, 47, 2, 45, 33, 47, 12, 12, 40, 61, 54, 9, 1, 30, 57, 42, 10, 14, 1, 22, 37, 25, 5, 36, 63, 24, 44, 34, 21, 56, 32, 18, 13, 59, 11]]


def rand_linear_permutation(n, density=0.5):
    """Returns the matrix representation of a linear permutation of
    {0,1}^n.

    This is done by generating random n x n binary matrices until one
    with full rank is found. As a random binary matrix has full rank
    with probability more than 1/2, this method is fine.

    """
    while True:
        result = [[0 for j in range(0, n)] for i in range(0, n)]
        for i, j in itertools.product(range(0, n), repeat=2):
            if random.random() < density:
                result[i][j] = 1
        result = Matrix(GF(2), n, n, result)
        if result.rank() == n:
            return result


def rand_linear_function(m, n, density=0.5):
    """Returns the matrix representation of a linear function mapping
    {0,1}^m to {0,1}^m.

    """
    result = [[0 for j in range(0, m)] for i in range(0, n)]
    for i, j in itertools.product(range(0, n), range(0, m)):
        if random.random() < density:
            result[i][j] = 1
    return Matrix(GF(2), n, m, result)


def linear_function_matrix_to_lut(mat):
    """Returns the codebook of the matrix mat."""
    result = [0 for x in range(0, 2**mat.ncols())]
    for x in range(0, 2**mat.ncols()):
        x_vec = to_bin(x, mat.ncols())
        y = from_bin(mat * vector(x_vec))
        result[x] = y
    return result


def print_bool(b):
    if b:
        print('-', end="")
    else:
        print('x', end="")


def test_ea_equivalence(s1, s2, nb_tests, test_from_code=True):
    """ Multiple tests between two EA affine mappings """
    n = int(log(len(s1), base=2))
    print("EA equivalence stress test. There should be no X, only -")
    # 1. The initial function are ea equivalent
    print_bool(are_ea_equivalent(s1, s2))
    if test_from_code:
        print_bool(are_ea_equivalent_from_code(s1, s2))
    # 2. Modify each of them by EA mappings and test again
    for _ in range(nb_tests):
        A1, A2, B1, B2 = [get_sbox(linear_function_matrix_to_lut(rand_linear_permutation(n))) for _ in range(4)]
        C1, C2 = [get_sbox(linear_function_matrix_to_lut(rand_linear_function(n, n))) for _ in range(2)]
        a1, a2, b1, b2 = [randint(0, len(s1)-1) for _ in range(4)]
        s1_bis = A1 * s1 * B1 + C1
        s2_bis = A2 * s2 * B2 + C2
        s1_bis = get_sbox([s1_bis[x ^ a1] ^ b1 for x in range(len(s1))])
        s2_bis = get_sbox([s2_bis[x ^ a2] ^ b2 for x in range(len(s2))])
        print_bool(are_ea_equivalent(s1_bis, s2_bis))
        if test_from_code:
            print_bool(are_ea_equivalent_from_code(s1_bis, s2_bis))
    print()


def test_random_ea_equivalence(n, nb_tests, test_from_code=True):
    print("Random function EA equivalence stress test. There should be no X, only -")
    for _ in range(nb_tests):
        f = [randint(0, 2**n-1) for _ in range(2**n)]
        g = [randint(0, 2 ** n - 1) for _ in range(2 ** n)]
        print_bool(not are_ea_equivalent(f, g))
        if test_from_code:
            print_bool(not are_ea_equivalent_from_code(f, g))
    print()


def test_banff_list():
    print("Banff list stress test. There should be no X, only -")
    for f in banff_list:
        for g in banff_list:
            if f != g:
                print_bool(not are_ea_equivalent(f, g))
                print_bool(not are_el_equivalent(f, g))
                print_bool(not are_ccz_equivalent(f, g))
                print_bool(not are_cczl_equivalent(f, g))
                print_bool(not are_affine_equivalent(f, g))
                print_bool(not are_linear_equivalent(f, g))
                print('|', end='')
        print()

def cardinal_of_automorphisms_group_banff_list():
    for f in banff_list:
        print(len(linear_equivalences(f, f)), end=" ")
        print(len(affine_equivalences(f, f)), end=" ")
        print(len(el_equivalences(f, f)), end=" ")
        print(len(ea_equivalences(f, f)), end=" ")
        print(len(cczl_equivalences(f, f)), end=" ")
        print(len(ccz_equivalences(f, f)), end=" ")
        print()


if __name__ == '__main__':
    n = 5
    ascon = get_sbox(list(sboxes.sboxes['Ascon']))

    # 1) THIS SHOULD WORK
    #test_ea_equivalence(ascon, ascon, 20)
    #test_random_ea_equivalence(n, 20)
    #test_banff_list()
    #cardinal_of_automorphisms_group_banff_list()

    #2) THERE IS STILL A REPRESENTATION PROBLEM BETWEEN THESE TWO FUNCTIONS
    # TODO Investigate
    aut = set([tuple(get_sbox(t).lut()) for t in automorphisms_from_ortho_derivative(banff_list[4])])
    aut_bis = set([tuple(get_sbox(t[0].transpose()).lut()) for t in ea_equivalences(banff_list[4], banff_list[4])])
    print(len(aut_bis), len(aut))
    print(aut == aut_bis)
