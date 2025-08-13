from sage.all import *
from sboxUv2 import *


n = 8
s = random_permutation_S_box(8)
print(differential_spectrum(s))
print(thickness_spectrum(s))
le_repr = le_class_representative(s)
print(le_repr)
print(differential_spectrum(le_repr))
