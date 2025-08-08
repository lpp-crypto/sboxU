from sage.all import *
from sboxUv2 import *

n = 4
g = GF(2**n)
cube = monomial(3, g)
print(cube)
print(differential_spectrum(cube))
print(ortho_derivative(cube))
