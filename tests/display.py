from sage.all import *
from sboxU import *

from sage.crypto.sboxes import sboxes

pprint(get_sbox(sboxes["AES"], name="AES"))

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
pi_prime = get_sbox(sboxes["Kuznyechik"])
pi_prime.rename("Kuznyechik")
lat_interactive_view(pi_prime)
pprint(fbct_spectrum(pi_prime))

F = get_sbox(sboxes["Skipjack"])
F.rename("Skipjack")
interactive_distribution_comparison_lat(F)
interactive_distribution_comparison_bct(F)
interactive_distribution_comparison_ddt(F, y_log_scale=False)
lat_interactive_view(F, with_sliders=False, with_cmap_choice=False)


S = random_permutation_S_box(8)
fbct_interactive_view(S)
pprint(fbct_spectrum(S))

fbct_interactive_view(monomial(3, GF(2**6)))
