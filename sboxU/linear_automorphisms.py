# Created on 2025-07-06 by baudrin-j.
# Modified on 2025-07-30 by baudrin-j.

from .sboxU_cython import *
from .utils import *
from .linear import *
from .ccz import table_linear_automorphisms
from time import time

# TODO add a test to compare with the linear automorphisms obtained from code equivalence.


def linear_automorphisms(lut, algorithm="alt_partition_diag_mappings", number_of_threads=8):
    """Return the linear automorphisms of a function F,
    i.e. all pairs of linear bijections (A, B) satisfying B o F o A = F.

    If the LAT is already computed and stored, the function linear_automorphisms_from_lat behaves the same
    but avoid the unnecessary recomputation of the LAT.

    More information on the search algorithms can be found in the corresponding C++ code.

    Args:
        lut: The look-up table of F.
        algorithm (str): A string describing the search algorithm to use. algorithm must be chosen among:
                    - "alt_partition_diag_mappings" (default)
                    - "alt_partition"
                    - "std_partition_diag_mappings"
        number_of_threads (int): The number of threads to use. By default, it is set to 8.

    Returns:
        list: A list of pairs of look-up tables corresponding to all pairs of linear bijections (A, B) satisfying B o F o A = F.

    """
    return linear_automorphisms_from_lat(lat(lut), algorithm, number_of_threads)


def is_linearly_self_equivalent(lut, algorithm="alt_partition_diag_mappings", number_of_threads=8):
    """ Checks whether a function F is linearly self-equivalent, i.e,
     if it exists a non-trivial pair of linear bijections (A, B) satisfying B o F o A = F.

     If the LAT is already computed and stored, the function is_linearly_self_equivalent_from_lat behaves the same
    but avoid the unnecessary recomputation of the LAT.

    More information on the search algorithms can be found in the corresponding C++ code.

    Args:
        lut: The look-up table of F.
        algorithm (str): A string describing the search algorithm to use. algorithm must be chosen among:
                    - "alt_partition_diag_mappings" (default)
                    - "alt_partition"
                    - "std_partition_diag_mappings"
        number_of_threads (int): The number of threads to use. By default, it is set to 8.

    Returns:
        tuple: A pair of look-up tables corresponding to a non-trivial (A, B) satisfying B o F o A = F, *if it exists*.
        Return False otherwise.

    """
    return is_linearly_self_equivalent_from_lat(lat(lut), algorithm, number_of_threads)


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
