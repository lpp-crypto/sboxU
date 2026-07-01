# automorphisms_from_ortho_derivative — mode benchmark

Checks that `automorphisms_from_ortho_derivative(f, mode="product")` — which enumerates the
full graph automorphism group via the semidirect product G1 ⋊ G2 (EL automorphisms times
derivative automorphisms) — produces the same count as the default `mode="standard"`,
that every map returned by `mode="product"` is a genuine graph automorphism of f,
and compares wall-clock times for both modes.

Tested on all quadratic APN functions in the 6-bit database.


## Preamble

```python
from time import time

DB_PATH = sixBitAPNs()
```


## Count and automorphism verification

A map L is a graph automorphism of f iff applying it to the graph of f gives the
graph of a function XOR-equivalent to f (same lut up to a constant shift).

```python
with APNFunctions(DB_PATH) as db:
    entries = db.query_functions({"degree": 2})
pprint("Loaded {} quadratic APN functions".format(len(entries)))

t_std = t_prod = 0.0
for entry in entries:
    f     = entry["sbox"]
    label = "id={}".format(entry["id"])
    t = time(); std  = automorphisms_from_ortho_derivative(f, mode="standard"); t_std  += time() - t
    t = time(); prod = automorphisms_from_ortho_derivative(f, mode="product");  t_prod += time() - t
    if len(std) != len(prod):
        fail("{}: standard={} automorphisms, product={}".format(label, len(std), len(prod)))
        continue
    for L in prod:
        g = ccz_equivalent_function(f, L)
        if len(g) == 0 or g.lut() != f.lut():
            fail("{}: map is not a graph automorphism".format(label))
            break
    else:
        success("{}: {} automorphisms, all verified".format(label, len(prod)))
pprint("Total time — standard: {:.3f}s  product: {:.3f}s".format(t_std, t_prod))
```
