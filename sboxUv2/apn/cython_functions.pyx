# -*- python -*-

from sboxUv2.core import Sb
from sboxUv2.config import MAX_N_THREADS


# !SECTION! Invariants

# !SUBSECTION! Invariants themselves

def ortho_derivative(q):
    """Returns the ortho-derivative of the function corresponding to the S_boxable object `q`, as defined e.g. in [TiT:CanCouPer22].

    Args:
        `q`: an S_boxable object corresponding to an APN function.

    Returns:
        An S_box instance containing the ortho-derivative of `q`. If the ortho-derivative of `q` is actually not defined (e.g. if it is not a quadratic APN), then returns an empty S_box.
    
    """
    sb = Sb(q)
    result = S_box(name="π_{".encode("UTF-8") + sb.name() + b"}")
    (<S_box>result).set_inner_sbox(
        cpp_ortho_derivative((<S_box>sb).cpp_sb[0])
    )
    return result


def sigma_multiplicities(s, k):
    # !TODO! docstring for sigma_multiplicities 
    sb = Sb(s)
    result = Spectrum(name="σ-mult".encode("UTF-8"))
    result.set_inner_sp(
        cpp_sigma_multiplicities((<S_box>sb).cpp_sb[0], k, MAX_N_THREADS)
    )
    return result


# !SUBSECTION! Aggregated invariants

def apn_ea_mugshot(s):
    sb = Sb(s)
    return cpp_apn_ea_mugshot((<S_box>sb).cpp_sb[0], MAX_N_THREADS)


# !SECTION!  CCZ-equivalence class exploration

def enumerate_ea_classes_apn_quadratic(
        s,
        n_threads=MAX_N_THREADS
):
    sb = Sb(s)
    result = []
    i = 0
    for new_s in cpp_enumerate_ea_classes_quadratic_apn(
        (<S_box>sb).cpp_sb[0],
        MAX_N_THREADS
    ):
        new_sb = S_box(name=b"CCZ-" + sb.name() + b"_" + str(i).encode("UTF-8"))
        new_sb.set_inner_sbox(<cpp_S_box>new_s)
        result.append(new_sb)
        i += 1
    return result


def ccz_equivalent_quadratic_function(
        s,
        n_threads=MAX_N_THREADS
):
    sb = Sb(s)
    result = S_box(name=b"deg2-CCZ-" + sb.name())
    result.set_inner_sbox(cpp_ccz_equivalent_quadratic_function(
        (<S_box>sb).cpp_sb[0],
        n_threads
    ))
    return result
