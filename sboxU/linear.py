#!/usr/bin/sage

from sage.all import *

import itertools
import random
from hashlib import sha256
from collections import defaultdict

from .sboxU_cython import *
from .utils import oplus
from .display import pretty_spectrum
from .diff_lin import lat_zeroes

DEFAULT_N_THREADS  = 16


# !SECTION! Utils for dealing with linear functions

def tobin(x, n):
    return [(x >> i) & 1 for i in reversed(range(0, n))]

def frombin(v):
    y = 0
    for i in range(0, len(v)):
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
    """Let $x = \\sum x_i 2^i$ and p = [p_0, ... p_{n-1}]. Returns
    $y = \\sum x_{p_i} 2^i$.

    """
    n = len(p)
    bin_x = [(x >> i) & 1 for i in range(0, n)]
    return sum(int(bin_x[p[i]] << i) for i in range(0, n))

def swap_halves(x, n):
    """If n=2k, swaps the first k bits of x with its last k bits."""
    if n % 2 == 1:
        raise Exception("Can't cut a word of {} bits in 2!".format(n))
    mask = sum(1 << i for i in range(0, n/2))
    l = x >> (n/2)
    r = x & mask
    return (r << (n/2)) | l


# !SUBSECTION! Fast linear mappings

class FastLinearMapping:
    """A convenient class to apply linear function on integers
    (interpreting them as elements of F_2^n) in way that is time
    efficient.

    """

    def __init__(self, L, lsb_first=False):
        if isinstance(L, sage.matrix.matrix0.Matrix): # case of a matrix
            self.inner_matrix = L
            if not lsb_first:
                self.masks = [
                    sum(int(L[i,j]) << (L.nrows()-i-1) for i in range(0, L.nrows()))
                    for j in reversed(range(0, L.ncols()))
                ]
            else:
                self.masks = [
                    sum(int(L[i,j]) << i for i in range(0, L.nrows()))
                    for j in range(0, L.ncols())
                ]
            self.input_size  = L.ncols()
            self.output_size = L.nrows()
        else:
            if "base_ring" in dir(L): # case of a polynomial
                gf = L.base_ring()
                self.masks = [L(gf.fetch_int(1 << shift)).integer_representation()
                              for shift in range(0, gf.degree())]
            if isinstance(L, list): # case of a list of masks
                self.masks = L[:]
        # setting input and output sizes
        self.inner_matrix = None
        self.input_size  = len(self.masks)
        self.output_size = 1
        masks_or = 0
        for m in self.masks:
            masks_or = masks_or | m
        cover = 1
        while (masks_or & cover) != masks_or:
            cover = (cover << 1) | 1
            self.output_size += 1

            
    def init_inner_matrix(self):
        self.output_size = self.input_size # this only works for square matrices
        mat = [[0 for j in range(0, self.input_size)]
               for i in range(0, self.output_size)]
        for j, m in enumerate(self.masks):
            for i in range(0, self.output_size):
                mat[i][self.input_size - j - 1] = (m >> (self.output_size - i - 1)) & 1
        self.inner_matrix = Matrix(GF(2),
                                   self.input_size,
                                   self.output_size,
                                   mat)

 

    def transpose(self):
        if self.inner_matrix == None:
            self.init_inner_matrix()
        return FastLinearMapping(self.inner_matrix.transpose())
    
    def inverse(self):
        if self.inner_matrix == None:
            self.init_inner_matrix()
        return FastLinearMapping(self.inner_matrix.inverse())


    # !TODO! add __rmult__ and __lmult__ functions to allow the
    # !composition of FastLinearMappings with other ones and with
    # !binary matrices


    def __call__(self, x):
        """Returns the result of applying L to the integer x, intepreting it
        as a binary vector.

        """
        result = 0
        for i in range(0, len(self.masks)):
            if (x >> i) & 1 == 1:
                result = oplus(result, self.masks[i])
        return result

    def __str__(self):
        return self.inner_matrix.str()
        

def block_FastLinearMapping(maps):
    """Return a FastLinearMapping obtained by combining those in the
    `maps` list of lists.

    For example, if `maps` is [[A,B], [C,D]] (where A,B,C and D are
    FastLinearMappings or objects that can be cast to one), then
    returns the FastLinearMapping corresponding to the block matrix

    [ A  B ]
    [ C  D ]

    meaning that it raises an Exception if the length of the lists
    don't match.

    Warning! only works when A, B, C and D are square matrices! This is
    because of the line marked with "<-": the computation of the
    output size of a FastLinearMapping is "wrong" if its image is
    equal over a smaller number of bits, in which case the processing
    of the blocks will give unexpected results.

    """
    # pre-processing
    for i in range(0, len(maps)):
        row = []
        for j in range(0, len(maps[0])):
            if not isinstance(maps[i][j], FastLinearMapping):
                maps[i][j] = FastLinearMapping(maps[i][j])
    # building all the masks
    masks = []
    for j in reversed(range(0, len(maps[0]))):
        for k in range(0, len(maps[0][j].masks)):
            v = 0
            shift = 0
            for i in reversed(range(0, len(maps))):
                v = v | (maps[i][j].masks[k] << shift)
                shift += maps[i][j].input_size # <---
            masks.append(v)
    return FastLinearMapping(masks)


    
# !SUBSECTION! Linear functions and their LUT


def linear_function_lut_to_matrix(l):
    """Turns the look up table `l` of a linear function into the
    corresponding binary matrix.

    """
    n = int(log(len(l), 2))
    result = []
    for i in range(0, n):
        line = [(int(l[1 << (n-1-j)]) >> (n-1-i)) & 1 for j in range(0, n)]
        result.append(line)
    return Matrix(GF(2), n, n, result)


def affine_function_lut_to_offset_and_matrix(a):
    """Turns the look up table `a` of an affine function into a list
    [c, L], where L is the binary matrix corresponding to its linear
    part, and where `c = a[0]`

    """
    n = int(log(len(a), 2))
    result = []
    l = [oplus(a[0], y) for y in a]
    for i in range(0, n):
        line = [(int(l[1 << (n-1-j)]) >> (n-1-i)) & 1 for j in range(0, n)]
        result.append(line)
    return [a[0], Matrix(GF(2), n, n, result)]


def linear_function_matrix_to_lut(mat):
    """Returns the codebook of the matrix mat."""
    result = [0 for x in range(0, 2**mat.ncols())]
    for x in range(0, 2**mat.ncols()):
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
        [(b >> (n - 1 - i)) & 1 for i in range(0, n)]
        for b in reversed(basis)]).transpose()
    check = linear_function_matrix_to_lut(result)
    if check[0:len(v)] == v:
        return result
    else:
        raise "no such matrix"


def F_2t_to_space(basis, n):
    """Returns the matrix corresponding to a permutation of (F_2)^n such
    that F_2^t (i.e. the set of integers < 2^t) is mapped to the space
    with the given basis using the function apply_bin_mat()

    """
    full_basis = complete_basis(basis, n)
    return Matrix(GF(2), n, n, [
        [(v >> (n-1-j)) & 1 for j in range(0, n)]
        for v in reversed(full_basis)
    ]).transpose()


def orthogonal_basis(B, n):
    """Returns a basis of the subspace of (F_2)^n that is orthogonal to
    all the vectors in B, the idea being that B is itself the basis of a
    subspace.

    """
    result = []
    r = 0
    v = 1
    while r < n-len(B):
        is_ortho = True
        for b_i in B:
            if scal_prod(b_i, v) != 0:
                is_ortho = False
                break
        if is_ortho:
            new_result = result + [v]
            new_r = rank_of_vector_set(new_result)
            if new_r > r :
                r = new_r
                result.append(v)
        v += 1
    return result


# !SUBSECTION! Vector/affine space extraction


def extract_bases(z,
                  dimension,
                  word_length,
                  n_threads=DEFAULT_N_THREADS,
                  number=b"fixed dimension"):
    """Returns a list containing the Gaussian Jacobi basis of each vector
    space of dimension `dimension` that is contained in the list `z` of
    integers intepreted as elements of $\\F_2^n$ where $n$ is equal to
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
      found, its bases will be returned and its subspaces will be
      ignored.

    """
    if dimension == 1:
        return [[x] for x in z]
    if len(z) < n_threads or 2**word_length < n_threads or 2**dimension < n_threads:
        n_threads = 1
    if number not in [b"all dimensions", b"fixed dimension", b"just one"]:
        raise Exception(b"Unknown value for parameter `number` in extract_bases:" + number)
    result = extract_bases_fast(z,
                               dimension,
                               word_length,
                               n_threads,
                               number)
    if number == b"all dimensions":
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
                         number=b"fixed dimension"):
    """Returns a list containing the Gaussian Jacobi basis of each affine
    space of dimension `dimension` that is contained in the list `z` of
    integers intepreted as elements of $\\F_2^n$ where $n$ is equal to
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
    if number not in [b"all dimensions", b"fixed dimension", b"just one"]:
        raise Exception(b"Unknown value for parameter `number` in extract_affine_bases:" + number)
    result = extract_affine_bases_fast(z,
                                     int(dimension),
                                     int(word_length),
                                     int(n_threads),
                                     number)
    if number == b"all dimensions":
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


def vector_spaces_bases_iterator(z, dimension, length):
    z.sort()
    sorted_z = []
    mask = 1
    current_msb_list = []
    for x in z:
        new_mask = False
        while (x & mask) != x:
            mask = (mask << 1) | 1
            new_mask = True
        if new_mask:
            if len(current_msb_list) > 0:
                sorted_z.append(current_msb_list[:])
            current_msb_list = []
        current_msb_list.append(x)
    if len(current_msb_list) > 0:
        sorted_z.append(current_msb_list)
    return vector_spaces_bases_iterator_rec(sorted_z, dimension, length)
        
    
def vector_spaces_bases_iterator_rec(sorted_z,
                                     dimension,
                                     word_length):
    if dimension == 1:
        result = []
        for z in sorted_z:
            for x in z:
                yield [x]
    elif len(sorted_z) >= dimension:
        tot_left = sum(len(z) for z in sorted_z)
        size_threshold = 2**(dimension-1) - 1
        for msb_index in range(0, len(sorted_z)-dimension+1):
            tot_left -= len(sorted_z[msb_index])
            if tot_left >= size_threshold:
                for v in sorted_z[msb_index]:
                    new_z = []
                    next_tot = 0
                    for i in range(msb_index+1, len(sorted_z)):
                        extracted = extract_vector(sorted_z[i], v)
                        if len(extracted) > 0:
                            new_z.append(extracted)
                            next_tot += len(extracted)
                    if next_tot >= size_threshold and len(new_z) >= dimension-1:
                        for rest_of_the_basis in vector_spaces_bases_iterator_rec(
                                new_z,
                                dimension-1,
                                word_length):
                            yield [v] + rest_of_the_basis
    
    


def vector_spaces_bases_iterator_old(z,
                                dimension,
                                word_length):
    """Returns an iterator going through the bases of all the vector
    spaces of dimension `dimension` that are included in `z`.

    The words in `z` are assumed to be of length `word_length`, and to
    be sorted in increasing order.

    """
    if dimension == 1:
        print("bleeeeh", z)
        return [[x] for x in z if x != 0]
    # sorting elements of z by their MSB
    z_by_msb = []
    current_msb_list = []
    mask = 1
    for x in z:
        new_msb = False
        while (x & mask) != x:
            new_msb = True
            mask = (mask << 1) | 1
        if new_msb:
            if len(current_msb_list) > 0:
                z_by_msb.append(current_msb_list)
            current_msb_list = [x]
        else:
            current_msb_list.append(x)
    if len(current_msb_list) > 0:
        z_by_msb.append(current_msb_list)
    if len(z_by_msb) >= dimension:
        # iterating the extraction
        for msb_index in range(0, len(z_by_msb)):
            for v_0 in z_by_msb[msb_index]:
                new_z = []
                for i in range(msb_index+1, len(z_by_msb)):
                    new_z += extract_vector(z_by_msb[i], v_0)
                if len(new_z) >= 2**(dimension-1) - 1:
                    for rest_of_the_basis in vector_spaces_bases_iterator(
                            new_z,
                            dimension-1,
                            word_length):
                        yield [v_0] + rest_of_the_basis
    

            

# !SUBSECTION!  Vector space bases and their properties



# def rank_deficit_of_vector_set_is_at_most(V, target):
#    """Returns whether c-r>=target where c is the number of elements in V
#    and r is the rank of the matrix obtained by "stacking" the n-bit
#    binary representation of the numbers in V.
#
#    """
#    return rank_deficit_of_vector_set_is_at_most_cpp(V, target)



def extract_basis(v, N):
    """Returns a subset of v such that the span of these elements is at
    least as big as v. In particular, if v is a vector space, it returns a
    base.

    """
    dim = rank_of_vector_set(v)
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

    The output *starts* with `basis`.

    """
    if rank_of_vector_set(basis) != len(basis):
        raise Exception("in complete_basis: the input must be independent! input={}".format(basis))
    r = len(basis)
    e_i = 1
    while r < N:
        new_basis = basis + [e_i]
        new_r = rank_of_vector_set(new_basis)
        if new_r > r:
            basis = new_basis[:]
            r = new_r
        e_i += 1
    return basis


def complete_basis_reversed(basis, N):
    """Returns a list of length N spanning the space 0..2^N-1 which
    contains the list of integers `basis`. Assumes that the elements
    of `basis` are linearly independent.

    The output *ends* with `basis`.
    """
    if rank_of_vector_set(basis) != len(basis):
        raise Exception("in complete_basis: the input must be independent! input={}".format(basis))
    r = len(basis)
    e_i = 1
    while r < N:
        new_basis = [e_i] + basis
        new_r = rank_of_vector_set(new_basis)
        if new_r > r:
            basis = new_basis[:]
            r = new_r
        e_i += 1
    return basis



def matrix_from_masks(basis, N):
    """Returns an NxN binary matrix M such that M*(1 << i) = basis[i] for
    all i < len(basis).

    """
    b = basis + [0]*(N - len(basis))
    return Matrix(GF(2), N, N, [
        [(b[i] >> j) & 1 for j in reversed(range(0, N))]
        for i in range(0, N)
    ]).transpose()


def get_generating_matrix(basis, N):
    """Returns an NxN binary matrix M such that M*(1 << i) = basis[i] for
    all i < len(basis) and such that M has full rank.

    """
    b = complete_basis(basis, N)
    return Matrix(GF(2), N, N, [
        [(b[i] >> j) & 1 for j in reversed(range(0, N))]
        for i in range(0, N)
    ]).transpose()


def linear_span(basis, with_zero=True):
    result = []
    if with_zero:
        result.append(0)
    visited = defaultdict(int)
    visited[0] = 1
    for i in range(1, 2**len(basis)):
        x = 0
        for j in range(0, len(basis)):
            if (i >> j) & 1 == 1:
                x = oplus(x, basis[j])
        if visited[x] != 1:
            result.append(x)
            visited[x] = 1
    return result


def bin_mat_to_int(m):
    """Turns a binary matrix into an integer via a simple bijection."""
    result = 0
    n_rows, n_cols  = len(m), len(m[0])
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            result = (result << 1) | m[i][j]
    return result


# !SUBSECTION! Easy interaction with finite fields

def mult_ff(x, y, F):
    return (F.fetch_int(x) * F.fetch_int(y)).integer_representation()


def div_ff(x, y, F):
    return (F.fetch_int(x) / F.fetch_int(y)).integer_representation()


def pow_ff(x, a, F):
    return (F.fetch_int(x)**a).integer_representation()



# !SECTION! Tests

def test_fast_multiplier(verbose=False):
    print("testing fast linear mappings")
    all_good = True
    n, m = 8, 4
    for index_L in range(0, 10):
        m += 1
        L = rand_linear_function(n, m)
        L_map = FastLinearMapping(L)
        if verbose:
            print("--- " + str(L_map.input_size()) + str(L_map.output_size()))
        for index_x in range(0, 8):
            x = randint(1, 2**n-1)
            y = apply_bin_mat(x, L)
            y_prime = L_map(x)
            if verbose:
                print(str(y) + " " + str(y_prime))
            if y != y_prime:
                all_good = False
                break
    if verbose:
        if all_good:
            print("[SUCCESS]")
        else:
            print("[FAIL]")
    return all_good


def test_vector_spaces_bases_iterator():
    from random import shuffle
    N = 10
    d = 4
    random_set = range(1, 2**N)
    shuffle(random_set)
    random_set = random_set[0:int(len(random_set)/3)]
    bases_0 = extract_bases(random_set, d, N)
    bases_1 = []
    for b in vector_spaces_bases_iterator(random_set, d, N):
        bases_1.append(b)
    print(bases_0)
    print("\n")
    print(bases_1)


if __name__ == '__main__':
    # print(test_fast_multiplier())
    test_vector_spaces_bases_iterator()
