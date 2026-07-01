# automorphisms_from_ortho_derivative — mode benchmark

Checks that `automorphisms_from_ortho_derivative(f, mode="product")` — which enumerates the
full graph automorphism group via the semidirect product G1 ⋊ G2 (EL automorphisms times
derivative automorphisms) — produces the same count as the default `mode="standard"`,
that every map returned by `mode="product"` is indeed a graph automorphism of f,
and compares wall-clock times for both modes.

Tested on all quadratic APN functions in the 6-bit database.


## Preamble

```python
import os
from time import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH    = os.path.join(SCRIPT_DIR, "../sboxU/sboxU/scripts/apnDB/apn6.db")
```


## Count and automorphism verification

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
    bad = next(
        (L for L in prod
         if len(g := ccz_equivalent_function(f, L)) == 0
         or g.lut() != f.lut()),
        None
    )
    if bad is None:
        success("{}: {} automorphisms, all verified".format(label, len(prod)))
    else:
        fail("{}: map {} is not a graph automorphism".format(label, bad))
pprint("Total time — standard: {:.3f}s  product: {:.3f}s".format(t_std, t_prod))
```
