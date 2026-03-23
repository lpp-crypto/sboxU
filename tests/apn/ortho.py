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
        field = GF(2**n)
        X = PolynomialRing(field, "X").gen()
        g = field.gen()

        cube = get_sbox(X**3)
        not_cube = get_sbox(X**3 + g**11*X**6 + g*X**9)

        A = rand_linear_permutation(n)
        B = rand_linear_permutation(n)
        C = rand_linear_function(n, n)
        print(A)
        print()
        print(B)
        other_func = get_sbox(B) * cube * get_sbox(A) + get_sbox(C)

        section("functions considered")
        
        for f in [cube, other_func, not_cube]:
            pprint(f)
            pprint(differential_spectrum(f))
            pprint(absolute_walsh_spectrum(f))
            pprint(degree_spectrum(f))
            pprint(thickness_spectrum(f))
            print()

        section("testing EA equivalence of (X^3) to (B o X^3 o A + C)")

        c = 0
        successes = 0
        for L in ea_mappings_from_ortho_derivative(cube, other_func):
            if len(xor_equivalence(cube, ccz_equivalent_function(other_func, L))) == 0:
                print("[FAIL]")
            else:
                successes += 1
            c += 1
        if c == 0:
            print("no mapping found")
        else:
            print("total: {} success out of {} ({:03.2f})%".format(
                  successes,
                  c,
                  100 * float(successes) / c
            ))

            
        section("testing EA equivalence of (X^3) to unequivalent function")

        c = 0
        successes = 0
        for L in ea_mappings_from_ortho_derivative(cube, not_cube):
            if len(xor_equivalence(cube, ccz_equivalent_function(other_func, L))) == 0:
                print("[FAIL]")
            else:
                successes += 1
            c += 1
        if c == 0:
            print("no mapping found")
        else:
            print("total: {} success out of {} ({:03.2f})%".format(
                  successes,
                  c,
                  100 * float(successes) / c
            ))
