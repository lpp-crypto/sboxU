from sage.all import *
from sage.crypto.sboxes import sboxes

from sboxUv2 import *


for t in range(0, 10):
    v = [randint(-4, 4) for i in range(0, 10)]
    sp = Spctr(v)
    v.sort()
    print("\n", v)
    print(sp.maximum(), sp.keys(), [(k,sp[k]) for k in sp.keys()])

for t in range(0, 10):
    s = random_function_S_box(randint(3, 4), randint(3, 4))
    print(s)
    print("-- DDT")
    dif = differential_spectrum(s)
    print(dif.maximum(), dif.keys(), dif)
    for row in ddt(s):
        print(row)
    print("-- LAT")
    wal = walsh_spectrum(s)
    print(wal.maximum(), wal.keys(), wal)
    l = lat(s)
    for row in l:
        print(row)
    if s == invert_lat(l):
        print("LAT inversion successfull")
    else:
        print("LAT inversion [FAILED]")

for name in ["Kuznyechik", "Fantomas"]:
    print("\n\n", name)
    s = Sb(sboxes[name])
    for tab in ["DDT", "LAT", "BCT"]:
        print(
            tab,
            table_anomaly(s, tab),
            table_negative_anomaly(sboxes[name], tab)
        )
    print("u = ", differential_uniformity(s))
    for u in range(6, 20, 2):
        if is_differential_uniformity_smaller_than(s, u):
            print("u <= ", u)
        else:
            print("u > ", u)

