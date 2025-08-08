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
    for row in lat(s):
        print(row)

for name in ["Kuznyechik", "Fantomas"]:
    print("\n\n", name)
    s = Sb(sboxes[name])
    for tab in ["DDT", "LAT"]:
        print(
            tab,
            table_anomaly(s, tab),
            table_negative_anomaly(sboxes[name], tab)
        )

