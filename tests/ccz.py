from sage.all import *
from sboxUv2 import *


n = 8
s = random_permutation_S_box(8)
print(differential_spectrum(s))
print(thickness_spectrum(s))
le_repr = le_class_representative(s)
print(le_repr)
print(differential_spectrum(le_repr))

n = 6
cube = monomial(3, GF(2**n))
print(thickness_spectrum(cube), cube.get_input_length(), cube.get_output_length())
s_s = enumerate_ea_classes(cube)
print("tot: ", len(s_s))
for s in s_s:
    print(s, thickness_spectrum(s))
