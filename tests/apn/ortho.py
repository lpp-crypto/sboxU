from sage.all import *
from sboxUv2 import *

import itertools 

def rand_linear_permutation(n):
    imgs = []
    r = 0
    while r < n:
        x = randint(1, 2**n-1)
        r_prime = rank_of_vector_set(imgs + [x])
        if r_prime > r:
            imgs.append(x)
            r = r_prime
    return Blm(imgs)

    
def rand_linear_function(m, n):
    imgs = []
    while len(imgs) < m:
        x = randint(1, 2**n-1)
        imgs.append(x)
    return Blm(imgs)


if __name__ == "__main__":
    with Experiment("testing ortho-derivative-based functions"):

        section("initialization")

        n = 6
        cube = monomial(3, GF(2**n))
        A = rand_linear_permutation(n)
        B = rand_linear_permutation(n)
        C = rand_linear_function(n, n)
        print(A)
        print()
        print(B)
        other_func = Sb(B) * cube * Sb(A) + Sb(C)

        section("functions considered")
        
        for f in [cube, other_func]:
            pprint(f)
            pprint(differential_spectrum(f))
            pprint(absolute_walsh_spectrum(f))
            pprint(degree_spectrum(f))
            pprint(thickness_spectrum(f))
            print()

        section("testing EA equivalence")

        c = 0
        for x in ea_mappings_from_ortho_derivative(cube, other_func):
            print(x)
            print()
            c += 1
        print("total: ", c)
