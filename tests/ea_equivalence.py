from sage.all import *
from sage.crypto import sboxes
from sboxUv2 import Sb, to_bin, from_bin, are_code_equivalent
import random
import itertools

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


if __name__ == '__main__':
    n = 5
    ascon = Sb(list(sboxes.sboxes['Ascon']))
    print(are_code_equivalent(ascon.lut(), (ascon ** -1).lut(), 'CCZ'))
    print(are_code_equivalent(ascon.lut(), (ascon ** -1).lut(), 'EA'))
    print('-'*30)
    for i in range(10):
        A = Sb(linear_function_matrix_to_lut(rand_linear_permutation(n)))
        B = Sb(linear_function_matrix_to_lut(rand_linear_permutation(n)))
        C = Sb(linear_function_matrix_to_lut(rand_linear_function(n, n)))
        print(are_code_equivalent((A * ascon * B + C).lut(), ascon.lut(), 'EA'))
    print('-'*30)
    for i in range(100):
        a = [randint(0, 2**n-1) for _ in range(2**n)]
        b = [randint(0, 2 ** n - 1) for _ in range(2 ** n)]
        print(are_code_equivalent(a, b, 'EA'))
