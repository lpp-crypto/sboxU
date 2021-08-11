"""Contains many 10-bit APN functions"""

from sage.all import GF, PolynomialRing

# global variables of the module
N = 10
F = GF(2**N, name="a")
g = F.gen()
POLY_RING = PolynomialRing(F, "X")
X = POLY_RING.gen()


def poly_to_lut(p):
    s = []
    for x_i in range(0, 2**N):
        y = (p(F.fetch_int(x_i))).integer_representation()
        s.append(y)
    return s


def all_quadratics():
    return [
        poly_to_lut(X**3),
        poly_to_lut(X**9),
    ]
