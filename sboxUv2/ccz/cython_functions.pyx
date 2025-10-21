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


cdef class WalshZeroesSpaces:
    def __init__(self):
        self.cpp_wzs = new cpp_WalshZeroesSpaces()
        self.mappings = []

        
    def image_by(self, L):
        Lm = Blm(L)
        result = WalshZeroesSpaces()
        (<WalshZeroesSpaces>result).cpp_wzs[0] = (<WalshZeroesSpaces>self).cpp_wzs[0].image_by(
            (<BinLinearMap>Lm).cpp_blm[0]
        )
        return result

    
    def init_mappings(self):
        self.cpp_wzs[0].init_mappings()
        for m in self.cpp_wzs[0].mappings:
            L = BinLinearMap()
            (<BinLinearMap>L).cpp_blm[0] = m
            self.mappings.append(L)
        

    def get_mappings(self):
        if len(self.mappings) == 0:
            self.init_mappings()
        return self.mappings
        

    def thickness_spectrum(self):
        result = Spectrum(name=b"thickness")
        result.set_inner_sp(
            (<WalshZeroesSpaces>self).cpp_wzs[0].thickness_spectrum()
        )
        return result
        


def get_WalshZeroesSpaces(s, n_threads=MAX_N_THREADS):
    sb = Sb(s)
    result = WalshZeroesSpaces()
    (<WalshZeroesSpaces>result).cpp_wzs = new cpp_WalshZeroesSpaces(
        (<S_box>sb).cpp_sb[0],
        n_threads
    )
    result.init_mappings()
    return result



# !SECTION! Exploring a CCZ-equivalence class


def ccz_equivalent_function(s, L):
    """Applies a linear permutation to the graph of a function, and returns the function whose graph is the result.
    Assumes that the linear permutation is admissible. If not, returns an empty S-box.

    Args:
        s: an S_box-able object
        L: a BinLinearMap-able object

    Returns:
        If `L` is indeed admissible for `s`, then returns an `S_box` instance corresponding to the function whose graph is the image of that of `s` by `L`. Otherwise, returns an empty `S_box`.
    
    """
    sb = Sb(s)
    basis = Blm(L)
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

    

def EA_mapping(A, B, C):
    """Assume that A and B are full-rank linear applications, and C another (maybe not-full rank) linear application. Given a function F, these can be used to define G as an extended affine equivalent function of F. This function helps expressing the relationship between F and G in terms of graph.
    
    Args:
        A: a BinLinearMap-able object corresponding to a full rank linear application.
        B: a BinLinearMap-able object corresponding to a full rank linear application.
        C: a BinLinearMap-able object.

    Returns:
        A `BinLinearMap` corresponding to the linear permutation that must be applied to the graph of a function F in order to obtain the graph of the function B*F*A + C.
    
    """
    result = BinLinearMap()
    Ablm = Blm(A)
    Bblm = Blm(B)
    Cblm = Blm(C)
    result.cpp_blm[0] = cpp_EA_mapping(
        (<BinLinearMap>Ablm).cpp_blm[0],
        (<BinLinearMap>Bblm).cpp_blm[0],
        (<BinLinearMap>Cblm).cpp_blm[0]
    )
    return result
    

# !SECTION! LAT-based equivalence tests

# !SUBSECTION! All linear automorphisms

def equivalences_from_lat(
        lat1, lat2,
        single_non_trivial_answer=False,
        equivalence_type="linear",
        n_threads=MAX_N_THREADS
):
    """ Computes all equivalence relations between two functions F, G, from their linear approximation tables (LAT).

        single_non_trivial_answer (bool) : Determines whether the search should stop when at the first automorphism that is encountered.
        number_of_threads (int): The number of threads to use. By default, it is set to 8.
        equivalence_type (string): Determines the type of the searched automorphisms.
        If "linear", search for all pairs of linear bijections (A, B) satisfying B o F o A = F.
        If "extended-search", search for all (A, B, C) linear with A, B bijective satisfying B o F o A + C = F.
        If "ccz-linear", search for all *linear* bijection A such that A(graph(F)) = graph(G).

        Returns:
            list: A list of pairs of BinLinearMaps corresponding to all pairs of linear bijections (A, B) satisfying B o F o A = F.

    """
    res = []
    for blocks in cpp_equivalences_from_lat(lat1, lat2, single_non_trivial_answer, n_threads, equivalence_type.encode('ascii')):
        abcd = []
        for b in blocks:
            new_blm = BinLinearMap()
            new_blm.cpp_blm[0] = b
            abcd.append(new_blm)
        res.append(abcd)
    return res

def linear_equivalences_from_lat(lat1, lat2, single_non_trivial_answer=False, n_threads=MAX_N_THREADS):
    return equivalences_from_lat(lat1, lat2, single_non_trivial_answer, "linear", n_threads)

def extended_linear_equivalences_from_lat(lat1, lat2, single_non_trivial_answer=False, n_threads=MAX_N_THREADS):
    return equivalences_from_lat(lat1, lat2, single_non_trivial_answer, "extended-linear", n_threads)

def ccz_linear_equivalences_from_lat(lat1, lat2, single_non_trivial_answer=False, n_threads=MAX_N_THREADS):
    return equivalences_from_lat(lat1, lat2, single_non_trivial_answer, "ccz-linear", n_threads)


def linear_equivalences(lut1, lut2, single_non_trivial_answer=False, n_threads=MAX_N_THREADS):
    return linear_equivalences_from_lat(lat(lut1), lat(lut2), single_non_trivial_answer, n_threads)

def extended_linear_equivalences(lut1, lut2, single_non_trivial_answer=False, n_threads=MAX_N_THREADS):
    return extended_linear_equivalences_from_lat(lat(lut1), lat(lut2), single_non_trivial_answer, n_threads)

def ccz_linear_equivalences(lut1, lut2, single_non_trivial_answer=False, n_threads=MAX_N_THREADS):
    return ccz_linear_equivalences_from_lat(lat(lut1), lat(lut2), single_non_trivial_answer, n_threads)


def are_linear_equivalent(lut1, lut2, n_threads=MAX_N_THREADS):
    t = linear_equivalences(lut1, lut2, True, n_threads)
    return lut1 == lut2 or len(t) > 0

def are_extended_linear_equivalent(lut1, lut2, n_threads=MAX_N_THREADS):
    t = extended_linear_equivalences(lut1, lut2, True, n_threads)
    return lut1 == lut2 or len(t) > 0

def are_ccz_linear_equivalent(lut1, lut2, n_threads=MAX_N_THREADS):
    t = ccz_linear_equivalences(lut1, lut2, True, n_threads)
    return lut1 == lut2 or len(t) > 0
