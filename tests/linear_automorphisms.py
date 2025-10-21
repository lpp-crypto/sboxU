from sage.all import *
from sboxUv2 import *
from time import time

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


# TODO Test consistency with the previous table_linear_automorphisms function of SboxUv1
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
        number_of_threads (int): The number of threads to use. By default, it is set to 8.
    """
    tab = lat(lut)
    automorphisms_algo_by_algo = []

    if algorithms is None:
        algorithms = ["alt_partition_diag_mappings", "alt_partition", "std_partition_diag_mappings"]

    print('Search time ', end='')
    for algo in algorithms:
        start = time()
        automorphisms = linear_automorphisms_from_lat(tab, algo, number_of_threads)
        automorphisms_algo_by_algo.append(automorphisms)
        stop = time()
        print("%.3f " % (stop - start), end='')
        Slut = Sb(lut)
        # Each output is indeed an automorphism
        for a, b in automorphisms:
            a = Sb(a.transpose().inverse())
            b = Sb(b.transpose())
            assert b * Slut * a == Slut
    print()

    # All algorithms return the same output
    for l in automorphisms_algo_by_algo[1:]:
        assert len(l) == len(automorphisms_algo_by_algo[0])
        for x in automorphisms_algo_by_algo[0]:
            assert x in l


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
    print("[alt_partition_diag_mappings, alt_partition, std_partition_diag_mappings]")

    # The C++ searches
    for algo in ["alt_partition_diag_mappings", "alt_partition", "std_partition_diag_mappings"]:
        start = time()
        automorphism = is_linearly_self_equivalent_from_lat(tab, algo, number_of_threads)
        stop = time()
        automorphism_algo_by_algo.append(automorphism)
        print("%.3f " % (stop - start), end='')

    if automorphism_algo_by_algo[0] is False:
        for ans in automorphism_algo_by_algo:
            assert ans is False
    else:
        Slut = Sb(lut)
        for ans in automorphism_algo_by_algo:
            # Each output is indeed an automorphism
            assert type(ans) is tuple
            a, b = ans
            a = Sb(a.transpose().inverse())
            b = Sb(b.transpose())
            id_lut = list(range(len(lut)))
            assert a.lut() != id_lut or b.lut() != id_lut
            assert b * Slut * a == Slut


if __name__ == "__main__":
    cubic_BL = [0, 0, 0, 1, 0, 2, 4, 7, 0, 4, 6, 3, 8, 14, 10, 13, 0, 8, 16, 25, 5, 15, 17, 26, 32, 44, 54, 59, 45, 35, 63, 48, 0, 16, 26, 36, 34, 48, 60, 0, 45, 57, 49, 11, 7, 17, 31, 39, 43, 28, 14, 23, 12, 57, 45, 54, 38, 21, 5, 24, 9, 56, 46, 49]
    for i_f, f in enumerate(all_quadratics_6() + [cubic_BL]):
        Slut = Sb(f)
        lat_f = lat(f)
        ccz_equivalences = ccz_linear_equivalences_from_lat(lat_f, lat_f, False, 8)
        el_equivalences = extended_linear_equivalences_from_lat(lat_f, lat_f, False, 8)
        lin_equivalences = linear_equivalences_from_lat(lat_f, lat_f, False, 8)
        print('Number of equivalences', len(ccz_equivalences), len(el_equivalences), len(lin_equivalences))
        for abcd in ccz_equivalences:
            # Transpose
            a, b, c, d = [Sb(x.transpose()) for x in abcd]
            c, d = d, c  # The top-right and bottom-left block are switched by the transpose
            # Rename
            a, b, c, d = a**(-1), b, c*a**(-1), d
            if d == [0 for _ in range(2**6)]:
                assert (b * Slut * a) + c == Slut
            else:
                print('ccz eq')
        for abcd in el_equivalences:
            assert abcd in cczgit_equivalences
        for abcd in lin_equivalences:
            assert abcd in lin_equivalences
        # Test of the function are_ccz_linear_equivalent on the known APN functions in dim 6
        for i_h, h in enumerate(all_quadratics_6() + [cubic_BL]):
            assert (f == h) == are_ccz_linear_equivalent(f, h)