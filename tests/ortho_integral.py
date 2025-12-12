from sage.all import *
from sboxUv2 import *

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


def func_to_vec(s):
    sb = Sb(s)
    result = []
    for x in sb.input_space():
        result += to_bin(sb[x], sb.get_output_length())
    return result


if __name__ == "__main__":
    with Experiment("implementing ortho-integration"):
        for q in ccz_class_representatives:
            o = ortho_derivative(q)
            E = F2LinearSystem(n*2**n)
            # ensuring pi(0) = 0 (thus removing constant addition)
            for i in range(0, n):
                E.add_equation([i])
            for i in range(0, n):
                L_vec = func_to_vec(g**i*X)
                E.remove_solution([x for x in range(0, n*2**n)
                                   if L_vec[x] == 1])
            for i in range(1, n):
                L_vec = func_to_vec(X**(2**i))
                E.remove_solution([x for x in range(0, n*2**n)
                                   if L_vec[x] == 1])
                
            for x, a in itertools.product(range(0, 2**n), range(1, 2**n)):
                if x != a:
                    eq = []
                    b = to_bin(o[a], n)
                    for i in range(0, n):
                        if b[i] == 1:
                            eq += [
                                i + n * x,          # F_i(x)
                                i + n * a,          # F_i(a)
                                i,                  # F_i(0)
                                i + n * oplus(x,a), # F_i(x+a)
                            ]
                    E.add_equation(eq)
            print(len(E.kernel()))
