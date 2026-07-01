# Product Walsh Match — Walsh Zero Space Orbit Test

Tests that the Walsh zero space orbit partition computed via `init_mappings(G1, G2)` logic
(using `product_walsh_match` as the inner check) yields the same number of orbits as
`enumerate_ea_classes_apn_quadratic`.  CCZ classes **2 and 3** of the 6-bit APN database.


## Preamble

```python
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH    = os.path.join(SCRIPT_DIR, "../sboxU/sboxU/scripts/apnDB/apn6.db")
CCZ_IDS    = [2, 3]

def linearize(L):
    return L + L.get_cste()

def span_of(vecs):
    s = {0}
    for b in vecs:
        s |= {x ^ b for x in s}
    return frozenset(s)

def basis_from_span(s):
    basis, covered = [], {0}
    for v in sorted(s):
        if v and v not in covered:
            basis.append(v)
            covered |= {x ^ v for x in covered}
    return basis

def make_uf(n):
    p = list(range(n))
    def find(x):
        while p[x] != x:
            p[x] = p[p[x]]; x = p[x]
        return x
    def union(a, b):
        a, b = find(a), find(b)
        if a != b: p[a] = b
    def partition():
        cls = {}
        for i in range(n):
            cls.setdefault(find(i), set()).add(i)
        return frozenset(frozenset(v) for v in cls.values())
    return find, union, partition
```


## Initialization

```python
with APNFunctions(DB_PATH) as db:
    seen_ccz = {cid: db.query_functions({"degree": 2, "ccz_id": cid})[0]["sbox"]
                for cid in CCZ_IDS}
pprint("Loaded CCZ classes: {}".format(sorted(seen_ccz.keys())))
```


## Two-group orbit partition

Implements `init_mappings(G1, G2)` in Python:

1. For every $g_2 \in G_2$ and every $V_i$, register $W = \mathrm{lin}(g_2)^{-T}(V_i)$ in
   the preimage map.  Two spaces that share a preimage are immediately merged.
2. For every preimage $W$, apply each $\mathrm{lin}(g_1)^T$ and merge the origins whenever
   the image lands in the preimage map.

`product_walsh_match(G1_t, G2_ti, W_basis, bases[j])` directly implements step 2 for a
specific target $V_j$; the loop below uses it to drive the union-find.

```python
two_part = {}
for cid in sorted(seen_ccz):
    f     = seen_ccz[cid]
    label = "CCZ-class {}".format(cid)
    wz    = get_WalshZeroesSpaces(f)
    bases = wz.get_bases()
    k     = len(bases)
    G1    = graph_el_automorphisms_from_ortho_derivative(f)
    G2    = graph_automorphisms_from_derivatives(f)
    G1_t  = [linearize(L).transpose() for L in G1]
    G2_ti = [linearize(L).transpose().inverse() for L in G2]
    find, union, partition = make_uf(k)
    # Step 1: preimage map — merge spaces that share a G2⁻ᵀ image
    preimage = {}
    for i, base in enumerate(bases):
        for g2 in G2_ti:
            W = span_of(g2(v) for v in base)
            if W in preimage:
                union(i, preimage[W])
            preimage[W] = find(i)
    # Step 2: for each preimage W, check if any G1ᵀ(W) is also a preimage
    for W_span, i in list(preimage.items()):
        W_basis = basis_from_span(W_span)
        for j, base_j in enumerate(bases):
            if find(i) != find(j):
                if len(product_walsh_match(G1_t, G2_ti, W_basis, base_j)) > 0:
                    union(i, j)
    two_part[cid] = partition()
    reps = sorted(min(orb) for orb in two_part[cid])
    success("{}: {} Walsh zero spaces → {} orbits, representatives: {}".format(
        label, k, len(two_part[cid]), reps))
```


## Comparison with enumerate_ea_classes_apn_quadratic

```python
for cid in sorted(seen_ccz):
    f        = seen_ccz[cid]
    label    = "CCZ-class {}".format(cid)
    ea       = enumerate_ea_classes_apn_quadratic(f)
    n_orbits = len(two_part[cid])
    n_ea     = len(ea)
    if n_orbits == n_ea:
        success("{}: {} orbits == {} EA classes".format(label, n_orbits, n_ea))
    else:
        fail("{}: {} orbits != {} EA classes".format(label, n_orbits, n_ea))
```
