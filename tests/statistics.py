from sage.all import *
from sboxUv2 import *


for t in range(0, 10):
    v = [randint(-4, 4) for i in range(0, 10)]
    sp = Spctr(v)
    v.sort()
    print("\n", v)
    print(sp.maximum(), sp.keys(), [(k,sp[k]) for k in sp.keys()])

    s = random_permutation_S_box(4)
    print(s)
    dif = differential_spectrum(s)
    print(dif.maximum(), dif.keys(), dif)
