# cython_functions.pyx
from ..sbox.cython_functions cimport S_box, cpp_S_box
from ...statistics.cython_functions cimport Spectrum, cpp_Spectrum
from .cython_functions cimport cpp_degree_spectrum,cpp_anf_component
from sage.all import GF,PolynomialRing,prod


def degree_spectrum(S_box s):
    py_result = Spectrum(name=b"Degree")  
    py_result.set_inner_sp(
        cpp_degree_spectrum(<const cpp_S_box&>s.cpp_sb[0])
    )
    return py_result

def algebraic_normal_form_coordinate(S_box s):
    R = PolynomialRing(GF(2), ['x{}'.format(i) for i in range(s.get_input_length())])
    coeffs = cpp_anf_component(<const cpp_S_box&>s.cpp_sb[0])
    vars = R.gens()
    P = R(0)
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        # obtenir la représentation binaire de l'indice
        bits = [(i >> j) & 1 for j in range(s.get_input_length())]
        monome = prod(var**b for var, b in zip(vars, bits))
        P += monome
        
    return P

def algebraic_normal_form(S_box s): #Returns a list of the anf of the coordinates
    return [algebraic_normal_form_coordinate(s.coordinate(i)) for i in range(s.get_output_length())]

