from sage.all import *
from sboxUv2 import *


n = 4
g = GF(2**n)
cube = monomial(3, g)
print(cube)
pprint(differential_spectrum(cube))
pprint(ortho_derivative(cube))
pprint(sigma_multiplicities(cube, 4))
pprint(thickness_spectrum(cube))

print(apn_ea_mugshot(cube))
