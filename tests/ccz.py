from sage.all import *
from sboxUv2 import *


n = 10
g = GF(2**n)
cube = monomial(3, g)
print(cube)
print(differential_spectrum(cube))
print(thickness_spectrum(cube))
