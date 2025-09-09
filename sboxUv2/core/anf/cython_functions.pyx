# -*- python -*-
# cython_functions.pyx

from ..sbox cimport *
from ..spectrum cimport *

from sboxUv2.core.f2functions import to_bin, from_bin
from sboxUv2.core.sbox import Sb


from sage.all import GF, PolynomialRing, prod


def degree_spectrum(s):
    """The degree spectrum describes the number of components of an S-box that have algebraic degree exactly k.

    Args:
       s: an S_box-able object

    Returns:
       A Spectrum object d such that d[k] is the number of components of `s` with algebraic degree exactly equal to k.
    """
    sb = Sb(s)
    py_result = Spectrum(name=b"Degree")  
    py_result.set_inner_sp(
        cpp_degree_spectrum((<S_box>s).cpp_sb[0])
    )
    return py_result


def algebraic_normal_form_coordinate(s, polynomial_vars=None):
    sb = Sb(s)
    if polynomial_vars == None:
        R = PolynomialRing(
            GF(2),
            ['x{}'.format(i) for i in range(sb.get_input_length())]
        )
        polynomial_vars = R.gens()
    elif (len(polynomial_vars) != s.get_input_length()):
        raise Exception("Wrong number of variables in `polynomial_vars`")
    coeffs = cpp_anf_component((<S_box>sb).cpp_sb[0])
    P = R(0)
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        monomial = prod(var**b for var, b in zip(polynomial_vars,
                                                 to_bin(i,sb.get_input_length())))
        P += monomial
    return P


def algebraic_normal_form(s, polynomial_vars=None):
    sb = Sb(s)
    return [algebraic_normal_form_coordinate(sb.coordinate(i),
                                             polynomial_vars=polynomial_vars)
            for i in range(sb.get_output_length())]


def eval_anf(anf, x):
    return from_bin(anf(to_bin(x)))


def eval_vect_anf(anfs, x):
    x_bin = to_bin(x)
    y = 0
    for i in range(0, len(anfs)):
        y = (<int>(anfs[i](x_bin)) << i) | y
    return y
