"""Contains many 7-bit APN functions"""

from sage.all import GF, PolynomialRing

# global variables of the module
N = 7
F = GF(2**N, name="a")
g = F.gen()
POLY_RING = PolynomialRing(F, "X")
X = POLY_RING.gen()


def poly_to_lut(p):
    s = []
    for x_i in xrange(0, 2**N):
        y = (p(F.fetch_int(x_i))).integer_representation()
        s.append(y)
    return s


def all_quadratics():
    """All the functions in Table 7 of
    http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.215.5432&rep=rep1&type=pdf
    
    (minus 2 which are not quadratic)
    """
    return [
        poly_to_lut(X**3),
        poly_to_lut(X**3 + sum(X**(9*2**i) for i in xrange(0, N))),
        poly_to_lut(X**34 + X**18 + X**5),
        poly_to_lut(X**3 + X**17 + X**33 + X**34),
        poly_to_lut(X**5),
        poly_to_lut(X**9),
        poly_to_lut(X**65 + X**10 + X**3),
        poly_to_lut(X**3 + X**9 + X**18 + X**66),
        poly_to_lut(X**3 + X**12 + X**17 + X**33),
        poly_to_lut(X**3 + X**17 + X**20 + X**34 + X**66),
    ]

def all_non_quadratics():
    return [
        poly_to_lut(X**13),
        poly_to_lut(X**57),
        poly_to_lut(X**(2**N-2)),
    ]
