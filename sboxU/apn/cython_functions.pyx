# -*- python -*-

from sboxU.cython_types cimport *
from sboxU.core import get_sbox
from sboxU.config import MAX_N_THREADS


from cython.operator cimport dereference


# !SECTION! Ortho-derivatives and friends


def ortho_derivative(q):
    """Returns the ortho-derivative of the function corresponding to the S_boxable object `q`, as defined e.g. in [TiT:CanCouPer22].

    Args:
        `q`: an S_boxable object corresponding to an APN function.

    Returns:
        An S_box instance containing the ortho-derivative of `q`. If the ortho-derivative of `q` is actually not defined (e.g. if it is not a quadratic APN), then returns an empty S_box.
    
    """
    sb = get_sbox(q)
    result = S_box(name="π_{".encode("UTF-8") + sb.name() + b"}")
    (<S_box>result).set_inner_sbox(
        cpp_ortho_derivative((dereference((<S_box>sb).cpp_sb)))
    )
    return result


def ortho_integral(s):
    """Returns the ortho-integral of the function corresponding to the S_boxable object `s`, as defined in a paper that is yet to be put anywhere (work in progress).

    Args:
        `s`: an S_boxable object.

    Returns:
        An S_box instance containing the ortho-integral of `s`, that is, to a quadratic APN function that has `s` as its ortho-derivative. If there is no such function, returns the empty S-box.
    
    """
    sb = get_sbox(s)
    result = S_box(name="∫_{".encode("UTF-8") + sb.name() + b"}")
    (<S_box>result).set_inner_sbox(
        cpp_ortho_integral(dereference((<S_box>sb).cpp_sb))
    )
    return result



# !SECTION! Invariants

# !SUBSECTION! Invariants themselves

def sigma_multiplicities(s, k=4):
    # !TODO! docstring for sigma_multiplicities 
    sb = get_sbox(s)
    result = Spectrum(name="σ-mult".encode("UTF-8"))
    result.set_inner_sp(
        cpp_sigma_multiplicities(dereference((<S_box>sb).cpp_sb), k, MAX_N_THREADS)
    )
    return result


# !SUBSECTION! Aggregated invariants

def apn_ea_mugshot(s):
    sb = get_sbox(s)
    return cpp_apn_ea_mugshot(dereference((<S_box>sb).cpp_sb), MAX_N_THREADS)


def apn_ea_mugshot_from_spectra(
        abs_walsh_spec,
        deg_spec,
        sig_mult,
        thk_spec
):
    return cpp_apn_ea_mugshot(
        dereference((<Spectrum>abs_walsh_spec).cpp_sp),
        dereference((<Spectrum>deg_spec).cpp_sp),
        dereference((<Spectrum>sig_mult).cpp_sp),
        dereference((<Spectrum>thk_spec).cpp_sp)
    )


# !SECTION! CCZ-equivalence class exploration


def automorphisms_from_ortho_derivative(s, n_threads=MAX_N_THREADS):
    """Returns all EA automorphisms of the quadratic APN function s.

    Each entry is a pair (L, delta) where L is the 2n×2n upper-triangular
    matrix encoding the automorphism and delta is the input shift.  See
    ccz_class.hpp for the full block-decomposition convention.

    To compare with the (r, a) pairs returned by ea_equivalences (which use
    a lower-triangular convention), use the canonical (A_affine, B) form::

        abcd_L = ccz_block_decomposition(L)
        A_affine_lut = [get_sbox(abcd_L[1]).inverse()[x] ^ delta
                        for x in range(2**n)]   # block_B^{-1} XOR delta
        B_lut        = get_sbox(abcd_L[0]).lut()  # block_A = L_B_T

        abcd_r = ccz_block_decomposition(r)
        A_affine_lut = [get_sbox(abcd_r[0])[x] ^ a
                        for x in range(2**n)]   # block_A XOR a
        B_lut        = get_sbox(abcd_r[1]).lut()  # block_B

    Args:
        s: an S-boxable object for a quadratic APN function.
        n_threads: number of threads. Defaults to MAX_N_THREADS.

    Returns:
        A list of (F2AffineMap, int) pairs (L, delta).
    """
    sb = get_sbox(s)
    result = []
    cdef std_vector[pair[cpp_F2AffineMap, BinWord]] automorphisms = \
        cpp_automorphisms_from_ortho_derivative(
            dereference((<S_box>sb).cpp_sb),
            n_threads
        )
    for L_and_delta in automorphisms:
        new_blm = F2AffineMap()
        (<F2AffineMap>new_blm).set_inner_map(L_and_delta.first)
        result.append((new_blm, int(L_and_delta.second)))
    return result


def ea_mappings_from_ortho_derivative(
        s,
        s_prime,
        n_threads=MAX_N_THREADS
):
    """Returns all the EL mappings L such that graph(s) = L(graph(s_prime)) + c, for some constant c that is not returned. Works only for quadratic APN functions since it is based on the ortho-derivative.

    Args:
        s: an S-boxable object corresponding to a quadratic APN function
        s_prime: an S-boxable object corresponding to a quadratic APN function
        n_threads: the number of threads to use. Defaults to `MAX_N_THREADS`.

    Returns:
        A list of F2AffineMaps L_i such that the graph of s is, up to a constant addition, the same as the image of the graph of s_prime under the linear permutation L_i.
    
    """
    sb, sb_prime = get_sbox(s), get_sbox(s_prime)
    result = []
    cdef std_vector[cpp_F2AffineMap] ea_mappings  = cpp_ea_mappings_from_ortho_derivative(
            dereference((<S_box>sb).cpp_sb),
            dereference((<S_box>sb_prime).cpp_sb),
            n_threads
    )
    for L in ea_mappings:
        new_blm = F2AffineMap()
        (<F2AffineMap>new_blm).set_inner_map(<cpp_F2AffineMap>L)
        result.append(new_blm)
    return result


def enumerate_ea_classes_apn_quadratic(
        s,
        n_threads=MAX_N_THREADS
):
    sb = get_sbox(s)
    result = []
    i = 0
    cdef std_vector[cpp_S_box] ea_classes = cpp_enumerate_ea_classes_quadratic_apn(
        dereference((<S_box>sb).cpp_sb),
        n_threads
    )
    for new_s in ea_classes:
        new_sb = S_box(name=b"CCZ-" + sb.name() + b"_" + str(i).encode("UTF-8"))
        new_sb.set_inner_sbox(<cpp_S_box>new_s)
        result.append(new_sb)
        i += 1
    return result


def ccz_equivalent_quadratic_function(
        s,
        n_threads=MAX_N_THREADS
):
    sb = get_sbox(s)
    result = S_box(name=b"deg2-CCZ-" + sb.name())
    result.set_inner_sbox(cpp_ccz_equivalent_quadratic_function(
        dereference((<S_box>sb).cpp_sb),
        n_threads
    ))
    return result

            
def get_WalshZeroesSpaces_quadratic_apn(s, n_threads=MAX_N_THREADS):
    sb = get_sbox(s)
    result = WalshZeroesSpaces()
    (<WalshZeroesSpaces>result).cpp_wzs = make_unique[cpp_WalshZeroesSpaces](
        dereference((<S_box>sb).cpp_sb),
        <unsigned int>n_threads
    )
    # handling the initilization of the mappings by hand
    dereference(result.cpp_wzs).init_mappings(
        cpp_automorphisms_from_ortho_derivative(dereference((<S_box>sb).cpp_sb),
                                                n_threads)
    )
    for m in dereference(result.cpp_wzs).mappings:
        L = F2AffineMap()
        (<F2AffineMap>L).set_inner_map(<cpp_F2AffineMap>m)
        result.mappings.append(L)
    return result



# !SECTION! Switching Neighbours

def non_trivial_sn(s,ne,ns):
    sb = get_sbox(s)
    result = []
    i = 0
    SW = cpp_non_trivial_sn(dereference((<S_box>sb).cpp_sb),<cpp_Integer> ne, <cpp_Integer> ns )
    for sw_u in SW: 
        res_u = []
        for new_s in sw_u:
            new_sb = S_box(name=b"SW-" + sb.name() + b"_" + str(i).encode("UTF-8"))
            new_sb.set_inner_sbox(<cpp_S_box>new_s)
            res_u.append(new_sb)
            i += 1
        result.append(res_u)
    return result
