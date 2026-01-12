from sage.all import *
from sboxUv2 import *
from time import *
field = GF(2**8)
X = PolynomialRing(field, "X").gen()
n=8

t = time()
plop = non_trivial_sn(0,Sb(X**3),4*2**n,4*(2**n)*(2**n))
print(time() - t)