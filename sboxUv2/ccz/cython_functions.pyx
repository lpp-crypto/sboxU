# -*- python -*-


from sboxUv2.core import Sb, Blm
from sboxUv2.config import MAX_N_THREADS
from sboxUv2.statistics import lat


# !SECTION! Walsh zeroes and friends


def thickness_spectrum(s, spaces=None):
    """Computes the thickness spectrum of the S_boxable object `s`.

    The thickness spectrum was introduced in [FFA:CanPer19] and is an extended-affine equivalence class invariant. Its computation requires the knowledge of the vector spaces of a specific dimension contained in the so-called Walsh zeroes of the function, meaning that this function relies on the vector space search algorithm from [AC:BonPerTia19].

    Args:
        - s: an S_boxable object.
        - spaces: if the Walsh zeroes of s are already known, then it is possible to pass the corresponding WalshZeroesSpaces object as a facultative argument to avoid recomputing it.

    Returns:
        A Spectrum instance such where the entry with key `k` is the number of spaces of dimension s.get_input_length() contained in the Walsh zeroes of `s` that have an intersection with V that is of dimension k.

    """
    sb = Sb(s)
    result = Spectrum(name=b"thickness")
    result.set_inner_sp(
        cpp_thickness_spectrum((<S_box>sb).cpp_sb[0], MAX_N_THREADS)
    )
    return result



# !SECTION! Exploring a CCZ-equivalence class

def ccz_equivalent_function(s, basis):
    sb = Sb(s)
    basis = Blm(basis)
    result = S_box(name=b"CCZ-" + sb.name())
    result.set_inner_sbox(
        cpp_ccz_equivalent_function(
            (<S_box>sb).cpp_sb[0],
            (<BinLinearMap>basis).cpp_blm[0]
            )
    )
    return result


def enumerate_ea_classes(s):
    sb = Sb(s)
    result = []
    i = 0
    for new_s in cpp_enumerate_ea_classes(
        (<S_box>sb).cpp_sb[0],
        MAX_N_THREADS):
        new_sb = S_box(name=b"CCZ-" + sb.name() + b"_" + str(i).encode("UTF-8"))
        new_sb.set_inner_sbox(<cpp_S_box>new_s)
        result.append(new_sb)
        i += 1
    return result


# !SECTION! LAT-based equivalence tests

# !SUBSECTION! All linear automorphisms

def linear_automorphisms_from_lat(lat, algorithm="alt_partition_diag_mappings", number_of_threads=8):
    """ Compute the linear automorphisms (*) of a function F from its linear approximation table (LAT).
        (*) i.e. all pairs of linear bijections (A, B) satisfying B o F o A = F.

        More information on the search algorithms can be found in the corresponding C++ code.

        Args:
            lat: The linear approximation table of F.
            algorithm (str): A string describing the search algorithm to use. algorithm must be chosen among:
                        - "alt_partition_diag_mappings" (default)
                        - "alt_partition"
                        - "std_partition_diag_mappings"
            number_of_threads (int): The number of threads to use. By default, it is set to 8.

        Returns:
            list: A list of pairs of look-up tables corresponding to all pairs of linear bijections (A, B) satisfying B o F o A = F.

    """
    return cpp_linear_automorphisms_from_lat(lat, algorithm.encode('ascii'), number_of_threads)


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


# !SUBSECTION! At most one linear automorphism

def is_linearly_self_equivalent_from_lat(lat, algorithm="alt_partition_diag_mappings", number_of_threads=8):
    """ Checks whether a function F is linearly self-equivalent (*) from its linear approximation table (LAT).
     (*) i.e. if it exists a non-trivial pair of linear bijections (A, B) satisfying B o F o A = F.

    More information on the search algorithms can be found in the corresponding C++ code.

    Args:
        lat: The linear approximation table of F.
        algorithm (str): A string describing the search algorithm to use. algorithm must be chosen among:
                    - "alt_partition_diag_mappings" (default)
                    - "alt_partition"
                    - "std_partition_diag_mappings"
        number_of_threads (int): The number of threads to use. By default, it is set to 8.

    Returns:
        A pair of look-up tables corresponding to a non-trivial (A, B) satisfying B o F o A = F, *if it exists*.
        Return False otherwise.

    """
    res = cpp_is_linearly_self_equivalent_from_lat(lat, algorithm.encode('ascii'), number_of_threads)
    if res == ([], []):
        return False
    else:
        return res

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