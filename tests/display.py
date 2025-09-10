from sage.all import *
from sboxUv2 import *

from sage.crypto.sboxes import sboxes

pprint(Sb([]))
pprint(Sb(sboxes["AES"], name="AES"))

f = random_function_S_box(5, 14)
pprint(f)
s = random_permutation_S_box(5)
print(s)
pprint(s)

d = differential_spectrum(s)
pprint(differential_spectrum(s))
pprint(walsh_spectrum(s))
pprint(absolute_walsh_spectrum(s))
pprint(thickness_spectrum(s))


print("Testing interactive tables display")
pi_prime = Sb(sboxes["Kuznyechik"])
pi_prime.rename("Kuznyechik")
lat_interactive_view(pi_prime)

F = Sb(sboxes["Skipjack"])
F.rename("Skipjack")
interactive_distribution_comparison_lat(F)



