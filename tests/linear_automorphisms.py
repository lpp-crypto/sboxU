from sage.all import *
from sboxUv2 import *


# global variables of the module
N = 6
F = GF(2**N, name="a")
g = F.gen()
POLY_RING = PolynomialRing(F, "X")
X = POLY_RING.gen()


def poly_to_lut(p):
    s = []
    for x_i in range(0, 2**N):
        if SAGE_VERSION < (9, 8):
            y = (p(F.fetch_int(x_i))).integer_representation()
        else:
            y = (p(F.from_integer(x_i))).to_integer()
        s.append(y)
    return s


def all_quadratics_6():
    """What follows is the Banff complete list of quadratic APN functions"""
    return [
        poly_to_lut(X**3),
        poly_to_lut(X**3 + g**11*X**6 + g*X**9),
        poly_to_lut(g*X**5 + X**9 + g**4*X**17 + g*X**18 + g**4*X**20 + g*X**24 + g**4*X**34 + g*X**40),
        poly_to_lut(g**7*X**3 + X**5 + g**3*X**9 + g**4*X**10 + X**17 + g**6*X**18),
        poly_to_lut(X**3 + g*X**24 + X**10), # 4 <-- KIM
        poly_to_lut(X**3 + g**17*(X**17 + X**18 + X**20 + X**24)),
        poly_to_lut(X**3 + g**11*X**5 + g**13*X**9 + X**17 + g**11*X**33 + X**48),
        poly_to_lut(g**25*X**5 + X**9 + g**38*X**12 + g**25*X**18 + g**25*X**36),
        poly_to_lut(g**40*X**5 + g**10*X**6 + g**62*X**20 + g**35*X**33 + g**15*X**34 + g**29*X**48),
        poly_to_lut(g**34*X**6 + g**52*X**9 + g**48*X**12 + g**6*X**20 + g**9*X**33 + g**23*X**34 + g**25*X**40),
        poly_to_lut(X**9 + g**4*(X**10 + X**18 ) + g**9*(X**12 + X**20 + X**40 )),
        poly_to_lut(g**52*X**3 + g**47*X**5 + g*X**6 + g**9*X**9 + g**44*X**12 + g**47*X**33 + g**10*X**34 + g**33*X**40),
        poly_to_lut(g*(X**6 + X**10 + X**24 + X**33) + X**9 + g**4*X**17),
    ]


def test_linear_automorphisms(lut, algorithms=None, number_of_threads=8):
    """ Test function for linear_automorphisms.
    Verifies that each result coincide for all algorithms and that the returned pairs indeed correspond to linear automorphisms.

    Tested for all algorithms, for number_of_threads= 1 and 8 on the ortho-derivatives of all quadratics functions from:
        - sboxU.known_functions.sixBitAPN
        - sboxU.known_functions.sevenBitAPN
        - some of the ones of sboxU.known_functions.eightBitAPN.

    Args:
        lut: The look-up table of F.
        algorithms (list): A list of strings describing the search algorithms to use. algorithm must be chosen among:
                    - "alt_partition_diag_mappings" (default)
                    - "alt_partition"
                    - "std_partition_diag_mappings"
                    - "sage"
        number_of_threads (int): The number of threads to use. By default, it is set to 8.
    """
    tab = lat(lut)
    automorphisms_algo_by_algo = []

    algo_sage = False
    if algorithms is None:
        algorithms = ["alt_partition_diag_mappings", "alt_partition", "std_partition_diag_mappings"]
        algo_sage = True
    else:
        if "sage" in algorithms:
            algo_sage = True
            algorithms.remove("sage")

    print(algorithms, end='')
    if algo_sage:
        print("sage")
    else:
        print()

    print('Search time ', end='')
    # The C++ searches
    for algo in algorithms:
        start = time()
        automorphisms = linear_automorphisms_from_lat(tab, algo, number_of_threads)
        stop = time()
        automorphisms = set([tuple([tuple(a), tuple(b)]) for a, b in automorphisms])
        automorphisms_algo_by_algo.append(automorphisms)
        print("%.3f " % (stop - start), end='')

    # The sage search
    if algo_sage:
        start = time()
        automorphisms = table_linear_automorphisms(tab)
        stop = time()
        automorphisms = set([tuple([tuple(a), tuple(b)]) for a, b in automorphisms])
        automorphisms_algo_by_algo.append(automorphisms)
        print("%.3f " % (stop - start), end='')

    print()

    print('Asserts time', end='')
    start = time()
    # All algorithms return the same output
    assert all([automorphisms == automorphisms_algo_by_algo[0] for automorphisms in automorphisms_algo_by_algo])
    # Each output is indeed an automorphism
    for a, b in automorphisms_algo_by_algo[0]:
        a = linear_function_lut_to_matrix(a).transpose().inverse()
        a = linear_function_matrix_to_lut(a)
        b = linear_function_lut_to_matrix(b).transpose()
        b = linear_function_matrix_to_lut(b)
        assert comp([b, lut, a]) == lut
    stop = time()
    print(" %.3f " % (stop - start))


def test_is_linearly_self_equivalent(lut, number_of_threads=8):
    """ Test function for is_linearly_self_equivalent.
    Verifies that the result provided by the C++ searches are in line with the sage search.

    Tested for all algorithms, for number_of_threads= 1 and 8 on the ortho-derivatives of all quadratics functions from:
        - sboxU.known_functions.sixBitAPN
        - sboxU.known_functions.sevenBitAPN
        - some of the ones of sboxU.known_functions.eightBitAPN

    Args:
        lut: The look-up table of F.
        number_of_threads (int): The number of threads to use. By default, it is set to 8.
    """
    tab = lat(lut)
    automorphism_algo_by_algo = []
    print("[alt_partition_diag_mappings, alt_partition, std_partition_diag_mappings, sage]")

    # The C++ searches
    for algo in ["alt_partition_diag_mappings", "alt_partition", "std_partition_diag_mappings"]:
        start = time()
        automorphism = is_linearly_self_equivalent_from_lat(tab, algo, number_of_threads)
        stop = time()
        automorphism_algo_by_algo.append(automorphism)
        print("%.3f " % (stop - start), end='')

    # The sage search
    start = time()
    automorphism_sage = table_linear_automorphisms(tab)
    stop = time()
    automorphism_sage = set([tuple([tuple(a), tuple(b)]) for a, b in automorphism_sage])
    print("%.3f " % (stop - start), end='')

    if len(automorphism_sage) == 1:
        assert all([automorphism is False for automorphism in automorphism_algo_by_algo])
    else:
        assert all([automorphism[0] or automorphism[1] for automorphism in automorphism_algo_by_algo])
    print()

if __name__ == "__main__":

    for i_f, funcs in enumerate([all_quadratics_6]):
        print()
        print('=== TESTING orthoderivatives of %d-bit functions ==='%(i_f+6))
        for f in funcs():
            #test_linear_automorphisms(ortho_derivative(f), algorithms=["alt_partition_diag_mappings"])
            test_is_linearly_self_equivalent(ortho_derivative(f))
