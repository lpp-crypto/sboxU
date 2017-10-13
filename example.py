#!/usr/bin/sage

from sage.all import *
from sboxU import random_permutation, affine_equivalence, oplus, rand_linear_permutation, apply_bin_mat

N = 5

# generating a random permutation g
g = random_permutation(N)

# generating the affine transformation
A = rand_linear_permutation(N)  # A is a SAGE Matrix instance
a = randint(0, 2**N-1)          # a is a random number of [0,2^n)
B = rand_linear_permutation(N)  # B is a SAGE Matrix instance
b = randint(0, 2**N-1)          # a is a random number of [0,2^n)
f = [oplus(apply_bin_mat(g[apply_bin_mat(oplus(x, a), A)], B), b)
     for x in xrange(0, 2**N)]
print "f = (B o g o A)(x + a) + b, where:"
print "A = \n", A.str()
print "a = ", a
print "B = \n", B.str()
print "b = ", b
print ""
# recovering the affine equivalence of f and g
ae = affine_equivalence(f, g)   # the algorithm by Biryukov et al
if len(ae) == 0:
    print "AE algorithm failed"
else:
    print "AE relation between f and g rediscovered!"
    print "f = (B o g o A)(x + a) + b, where:"
    print "A = \n", ae[0].str()
    print "a = ", ae[1]
    print "B = \n", ae[2].str()
    print "b = ", ae[3]
