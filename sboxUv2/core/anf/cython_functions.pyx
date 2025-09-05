# -*- python -*-
# cython_functions.pyx

from ..sbox cimport *
from ..spectrum cimport *

from ..sbox import Sb
from sage.all import GF, PolynomialRing, prod


def degree_spectrum(s):
    sb = Sb(s)
    py_result = Spectrum(name=b"Degree")  
    py_result.set_inner_sp(
        cpp_degree_spectrum((<S_box>s).cpp_sb[0])
    )
    return py_result

def algebraic_normal_form_coordinate(s):
    sb = Sb(s)
    R = PolynomialRing(
        GF(2),
        ['x{}'.format(i) for i in range(sb.get_input_length())]
    )
    coeffs = cpp_anf_component((<S_box>sb).cpp_sb[0])
    vars = R.gens()
    P = R(0)
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        # getting binary representation of the index
        # !QUESTION! why not use to_bin? 
        bits = [(i >> j) & 1 for j in range(sb.get_input_length())]
        monomial = prod(var**b for var, b in zip(vars, bits))
        P += monomial
    return P


def algebraic_normal_form(s):
    sb = Sb(s)
    return [algebraic_normal_form_coordinate(sb.coordinate(i))
            for i in range(sb.get_output_length())]

