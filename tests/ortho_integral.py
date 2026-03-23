from sage.all import *
from sboxU import *

import itertools
from collections import defaultdict


n = 6
field = GF(2**6)
X = PolynomialRing(field, "X").gen()
g = field.gen()

ccz_class_representatives = [
    # Banff list for quadratic functions
    Sb(X**3),
    Sb(X**3 + g**11*X**6 + g*X**9),
    Sb(g*X**5 + X**9 + g**4*X**17 + g*X**18 + g**4*X**20 + g*X**24 + g**4*X**34 + g*X**40),
    Sb(g**7*X**3 + X**5 + g**3*X**9 + g**4*X**10 + X**17 + g**6*X**18),
    Sb(X**3 + g*X**24 + X**10), # 4 <-- KIM
    Sb(X**3 + g**17*(X**17 + X**18 + X**20 + X**24)),
    Sb(X**3 + g**11*X**5 + g**13*X**9 + X**17 + g**11*X**33 + X**48),
    Sb(g**25*X**5 + X**9 + g**38*X**12 + g**25*X**18 + g**25*X**36),
    Sb(g**40*X**5 + g**10*X**6 + g**62*X**20 + g**35*X**33 + g**15*X**34 + g**29*X**48),
    Sb(g**34*X**6 + g**52*X**9 + g**48*X**12 + g**6*X**20 + g**9*X**33 + g**23*X**34 + g**25*X**40),
    Sb(X**9 + g**4*(X**10 + X**18 ) + g**9*(X**12 + X**20 + X**40 )),
    Sb(g**52*X**3 + g**47*X**5 + g*X**6 + g**9*X**9 + g**44*X**12 + g**47*X**33 + g**10*X**34 + g**33*X**40),
    Sb(g*(X**6 + X**10 + X**24 + X**33) + X**9 + g**4*X**17),
]


if __name__ == "__main__":
    with Experiment("Testing ortho-integration"):
        for q in ccz_class_representatives:
            o = ortho_derivative(q)
            s = ortho_integral(o)
            pprint(s)
            for f in [s, q]:
                pprint(differential_spectrum(f),
                       absolute_walsh_spectrum(f),
                       thickness_spectrum(f))
