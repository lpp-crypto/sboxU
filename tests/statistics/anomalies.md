# Anomalies


## Concept


### Testing

#### Setting up imports and parameters 

We obviously need to import `sboxU` along with `SAGE`. To simplify counting, we also import the `defaultdict` class from standard python module `collections`.

```{.python file=maximums.py}
from sage.all import *
from sboxUv2 import *
from collections import defaultdict
```

Statistical tests loose their relevance when the input size becomes too small so we fix `n=8`. Running a couple thousand tests is sufficient to gather enough information, hence `n_tested=2**13`.

```{.python file=maximums.py}
n = 8
n_tested = 2**13
```

#### LAT

```{.python file=maximums.py}
dis = defaultdict(int)
print("---- experimental")
for t in range(0, n_tested):
    dis[linearity(random_permutation_S_box(n))] += 1
for k in sorted(dis.keys()):
    if dis[k] != 0:
        print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))
print("---- theoretical")
dis = expected_linearity_distribution_permutation(n, n)
for k in sorted(dis.keys()):
    if dis[k] != 0.0:
        print("{:3d}: {:.4f}".format(k, dis[k]))
```
