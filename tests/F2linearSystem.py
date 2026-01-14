from sage.all import *
from sboxUv2 import *

import itertools
from collections import defaultdict


def basic_test():
    with Experiment("Testing the F2LinearSystem class"):
        E = F2LinearSystem(6)
        E.add_equation([0, 2])
        print(E)
        print("")
        E.add_equation([1, 3])
        print(E)
        print("")
        E.add_equation([0, 3])
        print(E)
        print("")
        # E.add_equation([0])
        # print(E)
        print(E.rank())
        E.add_equation([1,2,3, 4, 5])
        print(E)
        print(E.rank())
        print(E.kernel_as_bytes())    


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
    # # cubic function
    # Sb([0, 0, 0, 8, 0, 26, 40, 58, 0, 33, 10, 35, 12, 55, 46, 29, 0, 11, 12, 15, 4, 21, 32, 57, 20, 62, 18, 48, 28, 44, 50, 10, 0, 6, 18, 28, 10, 22, 48, 36, 8, 47, 16, 63, 14, 51, 62, 11, 5, 24, 27, 14, 11, 12, 61, 50, 25, 37, 13, 57, 27, 61, 39, 9])
]


            


def to_lut_coordinates(s):
    sb = Sb(s)
    result = []
    for x in sb.input_space():
        if sb[x] != 0:
            result.append(x)
    return result


# def from_vector(v, length):
#     result = []
#     for v_i in v:
#         for j in range(0, 8):
#             result.append((v_i >> j) & 1)
#     return result[:length]

        
if __name__ == "__main__":
    with Experiment("using switching neighbourgs to test linear solving"):
        n = 6
        gf = GF(2**n)

        if n > 6:
            ccz_class_representatives = [ monomial(3, GF(2**n)) ]

        for f in ccz_class_representatives[0:]:
            
            section("new function")
            pprint(differential_spectrum(ortho_derivative(f)))
            E = [
                F2LinearSystem(2**n)
                for u in range(0, 2**n)
            ]
    
            subsection("removing parasite solutions")
    
            for u in range(1, 2**n):
                # constant function
                E[u].remove_solution(to_lut_coordinates([1]*2**n))
                # linear functions
                for i in range(0, n):
                    E[u].remove_solution(to_lut_coordinates([
                        (x >> i) & 1 for x in range(0, 2**n)
                    ]))
                # # other coordinates
                L_u = generating_BinLinearMap([u], n)
                f_prime = Sb(L_u)**-1 * f
                for b in range(1, n):
                    E[u].remove_solution(to_lut_coordinates(
                        f_prime.coordinate(b) 
                    ))
    
                    
            subsection("building systems")
            for a in range(1, 2**n):
                for x,y in itertools.combinations(range(0, 2**n), r=2):
                    if (y != oplus(x, a)):
                        d_x = oplus(f[oplus(x, a)], f[x])
                        d_y = oplus(f[oplus(y, a)], f[y])
                        u = oplus(d_x, d_y)
                        E[u].add_equation([oplus(x, a), x, oplus(y, a), y])
    
            subsection("solving")
            counters = defaultdict(int)
            for u in range(1, 2**n):

                k = E[u].kernel_as_bits()
                counters[len(k)] += 1
                for g in k:
                    tot = []
                    for x in range(0, 2**n):
                        tot.append(oplus(u*g[x], f[x]))
                    # pprint(differential_spectrum(tot))
                    # print(algebraic_normal_form(g))
                    if not is_differential_uniformity_smaller_than(tot, 2):
                        print("[FAIL]")
                        pprint(differential_spectrum(tot))


            subsection("summary")
            pprint(counters)
        
