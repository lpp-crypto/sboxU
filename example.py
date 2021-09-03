#!/usr/bin/sage

from sage.all import *
from sboxU import *



def test_ea(N):
    # generating a random permutation g
    g = random_permutation(N)
    
    # generating the affine transformation
    A = rand_linear_permutation(N)  # A is a SAGE Matrix instance
    a = randint(0, 2**N-1)          # a is a random number of [0,2^n)
    B = rand_linear_permutation(N)  # B is a SAGE Matrix instance
    b = randint(0, 2**N-1)          # a is a random number of [0,2^n)
    f = [oplus(apply_bin_mat(g[apply_bin_mat(oplus(x, a), A)], B), b)
         for x in range(0, 2**N)]
    print("f = (B o g o A)(x + a) + b, where:")
    print("A = \n" + A.str())
    print("a = " + str(a))
    print("B = \n" + B.str())
    print("b = " + str(b))
    print("")
    # recovering the affine equivalence of f and g
    ae = affine_equivalence(f, g)   # the algorithm by Biryukov et al
    if len(ae) == 0:
        print("AE algorithm failed")
    else:
        print("AE relation between f and g rediscovered!")
        print("f = (B o g o A)(x + a) + b, where:")
        print("A = \n" + ae[0].str())
        print("a = " + str(ae[1]))
        print("B = \n" + ae[2].str())
        print("b = " + str(ae[3]))

        
def test_vector_extraction(N):
    z = []
    for x in range(0, 2**(N)+2**(N-1)):
        y = randint(0, 2**N-1)
        if y not in z:
            z.append(y)
    
    for number in ["just one", "all dimensions", "fixed dimension"]:
        #for number in ["all dimensions"]:
        print("\n\n=== {}\n".format(number))
        bases = extract_affine_bases(z, int(N/2), N, number=number)
        print("total = {}".format(len(bases)))
        for b in bases:
            print(str(len(b)) + "  " + pretty_vector(b))


if __name__ == '__main__':
    test_ea(6)
