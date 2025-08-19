# -*- python -*-


from sboxUv2.core import Sb
from sboxUv2.config import MAX_N_THREADS


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
    result = Spectrum()
    result.set_inner_sp(
        cpp_thickness_spectrum((<S_box>sb).cpp_sb[0], MAX_N_THREADS)
    )
    return result



# !SECTION! Exploring a CCZ-equivalence class

def ccz_equivalent_function(s, basis):
    sb = Sb(s)
    result = S_box(name=b"CCZ-" + sb.name())
    result.set_inner_sbox(
        cpp_ccz_equivalent_function(
            (<S_box>sb).cpp_sb[0],
            basis
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



