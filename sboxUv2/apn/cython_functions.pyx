# -*- python -*-

from sboxUv2.sbox import Sb
from sboxUv2.config import N_THREADS


def ortho_derivative(q):
    """Returns the ortho-derivative of the function corresponding to the S_boxable object `q`, as defined e.g. in [TiT:CanCouPer22].

    Args:
        - `q`: an S_boxable object corresponding to an APN function.

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
    sb = Sb(s)
    result = Spectrum()
    result.set_inner_sp(
        cpp_sigma_multiplicities((<S_box>sb).cpp_sb[0], k, N_THREADS)
    )
    return result
