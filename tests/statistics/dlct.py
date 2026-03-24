from sboxU.statistics import dlct, dlct_spectrum, dlct_uniformity
from sboxU import Spectrum
from sage.crypto.sboxes import sboxes

sb = sboxes["Midori_Sb0"]
n = 4
N = 1 << n

print("Testing DLCT structural properties on Midori_Sb0")

D = dlct(sb)

# For b=0: (-1)^{0 . v} = 1 for all v, so DLCT[a][0] = #{x} = N for all a.
if all(D[a][0] == N for a in range(N)):
    print("  DLCT[a][0] = N for all a: Success")
else:
    print("  DLCT[a][0] = N for all a: Failure")

# For a=0: F(x) XOR F(x XOR 0) = 0 for all x, so DLCT[0][b] = N for all b.
if all(D[0][b] == N for b in range(N)):
    print("  DLCT[0][b] = N for all b: Success")
else:
    print("  DLCT[0][b] = N for all b: Failure")

print("Testing dlct_spectrum on Midori_Sb0")

sp = dlct_spectrum(sb)

if isinstance(sp, Spectrum):
    print("  Returns Spectrum instance: Success")
else:
    print("  Returns Spectrum instance: Failure")

# The spectrum should account for all (N-1)^2 non-trivial (a,b) pairs.
n_pairs = sum(sp[k] for k in sp)
if n_pairs == (N - 1) ** 2:
    print("  Spectrum covers all non-trivial pairs: Success")
else:
    print("  Spectrum covers all non-trivial pairs: Failure ({} instead of {})".format(
        n_pairs, (N - 1) ** 2))

print("Testing dlct_uniformity on Midori_Sb0")

u = dlct_uniformity(sb)
u_manual = max(abs(D[a][b]) for a in range(1, N) for b in range(1, N))
if u == u_manual:
    print("  dlct_uniformity matches manual computation: Success")
else:
    print("  dlct_uniformity matches manual computation: Failure ({} vs {})".format(u, u_manual))
