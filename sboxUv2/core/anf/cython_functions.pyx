# -*- python -*-


from ..sbox cimport *
from ..spectrum cimport *

from sboxUv2.core.f2functions import to_bin, from_bin
from sboxUv2.core.sbox import Sb


from sage.all import GF, PolynomialRing, prod


# !SECTION! Information about the ANF


# !SUBSECTION! General information

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
        cpp_degree_spectrum((<S_box>sb).cpp_sb[0])
    )
    return py_result


def algebraic_degree(s):
    # !TODO! docstring for algebraic_degree
    sb = Sb(s)
    return cpp_algebraic_degree((<S_box>sb).cpp_sb[0])

# !SUBSECTION! Precise information on the degree

def is_degree_bigger_than(s,d):
    """
    Args :
        s : an S_box-able object
        d : the degree bound we want to test

    Returns:
        A boolean value, True if the algebraic degree of s is (strictly) bigger than d, False otherwise.
    """
    sb=Sb(s)
    return cpp_is_degree_bigger_than((<S_box>sb).cpp_sb[0],d)


# !SUBSECTION! The ANF itself


def algebraic_normal_form_coordinate(s, polynomial_vars=None):
    """The algebraic normal form of a boolean function is the multivariate polynomial representation of this function.

    Args:
       s: an S_box-able object of output length 1. 
       polynomial_vars: the variables used to express the polynomial. 

    Returns:
       A polynomial corresponding to the anf of the boolean function s, expressed in terms of polynomial_vars if not None.
    """
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
    P = polynomial_vars[0].parent()(0)
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        monomial = prod(var**b for var, b in zip(polynomial_vars,
                                                 to_bin(i,sb.get_input_length())))
        P += monomial
    return P


def algebraic_normal_form(s, polynomial_vars=None):
    """The algebraic normal form of a vectorial boolean function is the multivariate polynomial representation of the coordinates of this function.

    Args:
        s: an S_box-able object 
        polynomial_vars: the variables used to express the polynomial. 

    Returns:
        A list of polynomials corresponding to the anfs of the coordinate functions of s, expressed in terms of polynomial_vars if not None.
    """
    sb = Sb(s)
    if polynomial_vars == None:
        R = PolynomialRing(
            GF(2),
            ['x{}'.format(i) for i in range(sb.get_input_length())]
        )
        polynomial_vars = R.gens()
    elif (len(polynomial_vars) != s.get_input_length()):
        raise Exception("Wrong number of variables in `polynomial_vars`")
    return [algebraic_normal_form_coordinate(sb.coordinate(i),
                                             polynomial_vars=polynomial_vars)
            for i in range(sb.get_output_length())]


# !SECTION! Evaluating ANFs

def eval_anf(anf, x):
    """
    Evaluate a multivariate binary polynomial (in Algebraic Normal Form, ANF) 
    at a given integer input. The polynomial is defined on `n` variables. The integer `x` is interpreted 
    as an n-bit binary vector, where each bit corresponds to the value assigned 
    to one variable of the polynomial.

    Args:
        anf: a multivariate binary polynomial on n variables (in algebraic normal form). 
        x (int): an integer between 0 and 2**n - 1. Its n-bit binary representation is used as the variable assignment for the polynomial.

    Returns:
        The evaluation of the polynomial at the assignment given by `x` (always 0 or 1).
    """
    return anf(to_bin(x,len(anf.args())))


def eval_vect_anf(anfs, x):
    """
    Evaluate a vectorial Boolean function, represented as a list of multivariate 
    binary polynomials (in Algebraic Normal Form, ANF), at a given integer input. Each polynomial in `anfs` corresponds to one output coordinate of the function. 
    The integer `x` is interpreted as an n-bit binary vector, used as the variable 
    assignment for all polynomials. 

    Args:
        anfs (list): list of multivariate binary polynomials. All polynomials must have the same number of variables.  
        x (int): an integer whose n-bit binary representation defines the input assignment for the polynomials.

    Returns:
        the integer `y` representing the vector of polynomial evaluations, where bit i of `y` is the output of `anfs[i](x)`.
    """
    x_bin = to_bin(x,len(anfs[0].args()))
    y = 0
    for i in range(0, len(anfs)):
        y = (<int>(anfs[i](x_bin)) << i) | y
    return y
