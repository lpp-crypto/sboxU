SboxU v0.1


# Disclaimer

`sboxU` is at a very early stage. Its API might change over time!

# Description

SAGE/Python functions useful for studying S-boxes and Boolean
functions such as computing the DDT, computing the Walsh spectrum,
affine equivalence testing... Some of them are implemented in C++
and/or in parallel to improve their performances.


# Dependencies

The `SboxU` was only tested on Linux (Ubuntu 16.04). To install it,
you need the following packages.

```
libboost-python-dev
libpython-dev
sage
cmake
```

# Usage

Here is how to write a SAGE script using the functions from the
`sboxU` module. Put the `sboxU` folder in the same directory as your
script. Then, move to the `sboxU` folder and run:

```
cmake .
make
```

You can now import the `sboxU` python module (or any of its functions)
from your script in the parent folder of `sboxU`.


# Example
As an example of the functions provided by `sboxU`,
here is a script which tests the `affine_equivalence` function.

```
#!/usr/bin/sage

from sage.all import *
from sboxU import random_permutation, affine_equivalence, oplus, rand_linear_permutation, apply_bin_mat

N = 5
g = random_permutation(N)
f = random_permutation(N)
ae = affine_equivalence(f, g)
if len(ae) == 0:
    print "Random permutations are not affine equivalent"
else:
    print "Random permutations actually are affine equivalent!"
print ""

A = rand_linear_permutation(N)
a = randint(0, 2**N-1)
B = rand_linear_permutation(N)
b = randint(0, 2**N-1)
f = [oplus(apply_bin_mat(g[apply_bin_mat(oplus(x, a), A)], B), b)
     for x in xrange(0, 2**N)]
print "f = (B o g o A)(x + a) + b, where:"
print "A = \n", A.str()
print "a = ", a
print "B = \n", B.str()
print "b = ", b
print ""

ae = affine_equivalence(f, g)
if len(ae) == 0:
    print "AE algorithm failed"
else:
    print "AE relation between f and g found!"
    print "f = (B o g o A)(x + a) + b, where:"
    print "A = \n", ae[0].str()
    print "a = ", ae[1]
    print "B = \n", ae[2].str()
    print "b = ", ae[3]
```
