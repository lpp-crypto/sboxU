# -*- python -*-


from sboxUv2.core import Sb, get_F2AffineMap
from sboxUv2.config import MAX_N_THREADS
from sboxUv2.statistics import lat

from cython.operator cimport dereference



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
        cpp_thickness_spectrum(dereference((<S_box>sb).cpp_sb), MAX_N_THREADS)
    )
    return result


# !SUBSECTION! The  WalshZeroesSpaces class

cdef class WalshZeroesSpaces:
    def __init__(self):
        self.cpp_wzs = new cpp_WalshZeroesSpaces()
        self.mappings = []

        
    def __dealloc__(self):
        self.cpp_wzs[0].destruct()
        free(self.cpp_wzs)

        
    def image_by(self, L):
        Lm = get_F2AffineMap(L)
        result = WalshZeroesSpaces()
        (<WalshZeroesSpaces>result).cpp_wzs[0] = (<WalshZeroesSpaces>self).cpp_wzs[0].image_by(
            dereference((<F2AffineMap>Lm).cpp_map)
        )
        return result

    
    def init_mappings(self):
        self.cpp_wzs[0].init_mappings()
        for m in self.cpp_wzs[0].mappings:
            L = F2AffineMap()
            L.set_inner_map(m)
            self.mappings.append(L)
        
    
    # def init_mappings_using_automorphisms(
    #         self,
    #         std_vector[F2AffineMap] automorphisms
    # ):
    #     std_vector(
    #     self.cpp_wzs[0].init_mappings()
    #     for m in self.cpp_wzs[0].mappings:
    #         L = F2AffineMap()
    #         (<F2AffineMap>L).cpp_blm[0] = m
    #         self.mappings.append(L)
        

    def get_mappings(self):
        if len(self.mappings) == 0:
            self.init_mappings()
        return self.mappings
        

    def get_mappings(self):
        if len(self.mappings) == 0:
            self.init_mappings()
        return self.mappings

    
    def get_bases(self) -> list:
        return [b.get_basis() for b in self.cpp_wzs[0].bases]
            
        
    def thickness_spectrum(self):
        result = Spectrum(name=b"thickness")
        result.set_inner_sp(
            (<WalshZeroesSpaces>self).cpp_wzs[0].thickness_spectrum()
        )
        return result

    
    def __iter__(self):
        for b in self.cpp_wzs[0].bases:
            yield BinLinearBasis(b.get_basis())
        return
        

def get_WalshZeroesSpaces(s, n_threads=MAX_N_THREADS):
    sb = Sb(s)
    result = WalshZeroesSpaces()
    (<WalshZeroesSpaces>result).cpp_wzs = new cpp_WalshZeroesSpaces(
        dereference((<S_box>sb).cpp_sb),
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
        L: a F2AffineMap-able object

    Returns:
        If `L` is indeed admissible for `s`, then returns an `S_box` instance corresponding to the function whose graph is the image of that of `s` by `L`. Otherwise, returns an empty `S_box`.
    
    """
    sb = Sb(s)
    basis = get_F2AffineMap(L)
    result = S_box(name=b"CCZ-" + sb.name())
    result.set_inner_sbox(
        cpp_ccz_equivalent_function(
            dereference((<S_box>sb).cpp_sb),
            dereference(((<F2AffineMap>basis).cpp_map))
            )
    )
    return result


def enumerate_ea_classes(s):
    sb = Sb(s)
    result = []
    i = 0
    cdef std_vector[cpp_S_box] ea_classes = cpp_enumerate_ea_classes(
        dereference((<S_box>sb).cpp_sb), 
        MAX_N_THREADS
    )
    for new_s in ea_classes :
        new_sb = S_box(name=b"CCZ-" + sb.name() + b"_" + str(i).encode("UTF-8"))
        new_sb.set_inner_sbox(<cpp_S_box>new_s)
        result.append(new_sb)
        i += 1
    return result


def enumerate_permutations_in_ccz_class(s):
    sb = Sb(s)
    result = []
    i = 0
    cdef std_vector[cpp_S_box] permuations_in_ccz_class = cpp_enumerate_permutations_in_ccz_class(
        dereference((<S_box>sb).cpp_sb), 
        MAX_N_THREADS
    )
    for new_s in permuations_in_ccz_class :
        new_sb = S_box(name=b"CCZ-" + sb.name() + b"_" + str(i).encode("UTF-8"))
        new_sb.set_inner_sbox(<cpp_S_box>new_s)
        result.append(new_sb)
        i += 1
    return result

    

def EA_mapping(A, B, C):
    """Assume that A and B are full-rank linear applications, and C another (maybe not-full rank) linear application. Given a function F, these can be used to define G as an extended affine equivalent function of F. This function helps expressing the relationship between F and G in terms of graph.
    
    Args:
        A: a F2AffineMap-able object corresponding to a full rank linear application.
        B: a F2AffineMap-able object corresponding to a full rank linear application.
        C: a F2AffineMap-able object.

    Returns:
        A `F2AffineMap` corresponding to the linear permutation that must be applied to the graph of a function F in order to obtain the graph of the function B*F*A + C.
    
    """
    result = F2AffineMap()
    Ablm = get_F2AffineMap(A)
    Bblm = get_F2AffineMap(B)
    Cblm = get_F2AffineMap(C)
    result.set_inner_map(cpp_EA_mapping(
        dereference((<F2AffineMap>Ablm).cpp_map),
        dereference((<F2AffineMap>Bblm).cpp_map),
        dereference((<F2AffineMap>Cblm).cpp_map),
    ))
    return result
    

# !SECTION! LAT-based equivalence tests

# !SUBSECTION! All "linear" equivalences

def equivalences_from_lat(
        sbox1, sbox2,
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
            list: A list of pairs of F2AffineMaps corresponding to all pairs of linear bijections (A, B) satisfying B o F o A = F.

    """
    res = []
    cdef std_vector[cpp_F2AffineMap] equivalences_from_lat = cpp_equivalences_from_lat(
            dereference((<S_box>sbox1).cpp_sb), 
            dereference((<S_box>sbox2).cpp_sb), 
            single_non_trivial_answer,
            n_threads,
            equivalence_type.encode('ascii')
    )
    for solution in equivalences_from_lat :
        new_blm = F2AffineMap()
        (<F2AffineMap>new_blm).set_inner_map(solution)
        res.append(new_blm)
    return res

# !TODO! in all functions below:
# ! 1. write docstrings
# ! 2. wrap them in simpler tests with Boolean results

def linear_equivalences(sbox1, sbox2, single_non_trivial_answer=False, n_threads=MAX_N_THREADS):
    return equivalences_from_lat(
        Sb(sbox1),
        Sb(sbox2),
        single_non_trivial_answer,
        "linear",
        n_threads
    )


def el_equivalences(sbox1, sbox2, single_non_trivial_answer=False, n_threads=MAX_N_THREADS):
    return equivalences_from_lat(
        Sb(sbox1),
        Sb(sbox2),
        single_non_trivial_answer,
        "extended-linear",
        n_threads
    )


def cczl_equivalences(sbox1, sbox2, single_non_trivial_answer=False, n_threads=MAX_N_THREADS):
    return equivalences_from_lat(
        Sb(sbox1),
        Sb(sbox2),
        single_non_trivial_answer,
        "ccz-linear",
        n_threads
    )


# !SUBSECTION! All "affine" equivalences

def up_to_constant_equivalences(sbox1, sbox2, single_non_trivial_answer=False,  equivalence_type="linear", n_threads=MAX_N_THREADS):
    shifted_sbox1 = Sb([x ^ sbox1[0] for x in sbox1])

    results = []
    for a in range(len(sbox1)):
        shifted_sbox2 = Sb([sbox2[x ^ a] ^ sbox2[a] for x in range(len(sbox2))])
        res = equivalences_from_lat(shifted_sbox1,
                                    shifted_sbox2,
                                    False,
                                    equivalence_type,
                                    8)
        for r in res:
            verify_up_to_constant_equivalence(sbox1,
                                              sbox2,
                                              r,
                                              a,
                                              equivalence_type) # Sanity checks
            results.append((r, a))
    return results


def affine_equivalences(sbox1, sbox2, single_non_trivial_answer=False, n_threads=MAX_N_THREADS):
    return up_to_constant_equivalences(
        Sb(sbox1),
        Sb(sbox2),
        single_non_trivial_answer,
        "linear",
        n_threads
    )


def ea_equivalences(sbox1, sbox2, single_non_trivial_answer=False, n_threads=MAX_N_THREADS):
    return up_to_constant_equivalences(
        Sb(sbox1),
        Sb(sbox2),
        single_non_trivial_answer,
        "extended-linear",
        n_threads
    )


def ccz_equivalences(sbox1, sbox2, single_non_trivial_answer=False, n_threads=MAX_N_THREADS):
    return up_to_constant_equivalences(
        Sb(sbox1),
        Sb(sbox2),
        single_non_trivial_answer,
        "ccz-linear",
        n_threads
    )


# !SUBSECTION! Boolean tests

# !TODO! in all functions below:
# ! 1. write docstrings
# ! 2. handle invariant computations before launching the full search
# ! 2. handle option to use either equivalence from lats or from codes


def are_linear_equivalent(sbox1, sbox2, n_threads=MAX_N_THREADS):
    t = linear_equivalences(sbox1, sbox2, True, n_threads)
    return sbox1 == sbox2 or len(t) > 0

def are_el_equivalent(sbox1, sbox2, n_threads=MAX_N_THREADS):
    t = el_equivalences(sbox1, sbox2, True, n_threads)
    return sbox1 == sbox2 or len(t) > 0

def are_cczl_equivalent(sbox1, sbox2, n_threads=MAX_N_THREADS):
    t = cczl_equivalences(sbox1, sbox2, True, n_threads)
    return sbox1 == sbox2 or len(t) > 0

def are_affine_equivalent(sbox1, sbox2, n_threads=MAX_N_THREADS):
    t = affine_equivalences(sbox1, sbox2, True, n_threads)
    return sbox1 == sbox2 or len(t) > 0

def are_ea_equivalent(sbox1, sbox2, n_threads=MAX_N_THREADS):
    t = ea_equivalences(sbox1, sbox2, True, n_threads)
    return sbox1 == sbox2 or len(t) > 0

def are_ccz_equivalent(sbox1, sbox2, n_threads=MAX_N_THREADS):
    t = ccz_equivalences(sbox1, sbox2, True, n_threads)
    return sbox1 == sbox2 or len(t) > 0


# !SUBSECTION! Some utils

def ccz_block_decomposition(blm):
    """ Decompose a 2n x 2n matrix (F2AffineMap) into the 4 corresponding n x n matrices (returned as F2AffineMap)."""
    abcd = []
    cdef std_vector[cpp_F2AffineMap] blocks = cpp_ccz_block_decomposition(dereference((<F2AffineMap>blm).cpp_map))
    for block in blocks :
        new_blm = F2AffineMap()
        new_blm.set_inner_map(block)
        abcd.append(new_blm)
    return abcd


def verify_up_to_constant_equivalence(sbox1, sbox2, L, a, equivalence_type):
    """ Verify whether sbox1 and sbox2 are indeed equivalent given an affine mapping x \mapsto Lx + a. """
    if equivalence_type in ["extended-linear", "linear"]:
        A, B, C, D = ccz_block_decomposition(L)
        A, B, C, D = Sb(A), Sb(B), Sb(C), Sb(D)

        # Verify the null blocks
        assert D == [0 for _ in range(len(D))]
        if equivalence_type == "linear":
            assert C == [0 for _ in range(len(C))]

        A_aff = Sb([A[x] ^ a for x in range(len(A))])
        B_aff = Sb([B[x] ^ sbox2[a] ^ B[sbox1[0]] for x in range(len(B))])
        
        sbox1_prime = (B_aff * sbox1) + C
        sbox2_prime = sbox2 * A_aff
        assert sbox1_prime == sbox2_prime
    if equivalence_type == "ccz-linear":
        pass
        # !TODO! Add appropriate test for ccz equivalence
