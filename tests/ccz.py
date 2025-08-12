from sage.all import *
from sboxUv2 import *


n = 7
g = GF(2**n)
cube = monomial(3, g)
print(cube)
print(differential_spectrum(cube))
print(thickness_spectrum(cube))
le_repr = le_class_representative(cube)
print(le_repr)
print(differential_spectrum(le_repr))
