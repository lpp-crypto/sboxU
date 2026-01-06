# Anomalies


## Concept

An "anomaly" is a concept that was introduced in [AC:BonPerTia19] (though its principle was already present in [C:BirPer15]). It can be positive or negative, and it quantifies how good the properties of a given S-box are compared to a permutation picked uniformly at random from the relevant set. If a an S-box has a high positive anomaly for a given property, then it means that its property is much better than should be expected from a random S-box. If instead its negative anomaly is high, then its property is much worse than should be expected.


Formally, consider a property $P$ defined for any S-box $S$, such that its outputs can be strictly ordered (for instance, if it's a number then we can simply compare two outputs).



## Testing

Let's see how this notion behaves. To this end, we obviously need to import `sboxU` along with `SAGE`. To simplify counting, we also import the `defaultdict` class from standard python module `collections`.

```python
from sage.all import *
from sboxUv2 import *
from collections import defaultdict
```


### Computing some anomalies

We will look at some S-boxes from the literature. To this end, we grab the `sboxes` dictionary shipped with SAGE which contains many of them.

```python
from sage.crypto.sboxes import sboxes
```

We are now ready to make some experiments.

```python
with Experiment("Checking some anomalies"):
```

#### Skipjack

Let's start with the S-box of Skipjack and let's look at the anomalies corresponding to the main tables:

```python
    section("Skipjack")
    skipjack = Sb(sboxes["Skipjack"])
    for table in ["LAT", "BCT", "DDT"]:
        print("{} {:10.5f} {:10.5f}".format(
            table,
            table_anomaly(skipjack, table),
            table_negative_anomaly(skipjack, table)
        ))
```

As we can see, the positive anomaly for the LAT is very high (remember that it is a logarithm): the S-box has been somehow optimized to have good linear property (see [C:BirPer15] and [PhD:Perrin17] for more details). Another way this property can be seen is through the display of a graph of its absolute Walsh spectrum along with what is expected for a random permutation. This is only done using the following function, which will open a new window displaying this graph. The blue area also shows the discrepancy between the expected behaviour of a random S-box and that of Skipjack.

(you must close the graph window for the program to continue)

```python
    interactive_distribution_comparison_lat(skipjack)
```

#### Other S-boxes

We can look at these values for other S-boxes: those of the AES, Fantomas, and Kuznyechik.

```python
    section("Other S-boxes")
    for k in ["AES", "Fantomas", "Kuznyechik"]:
        s = Sb(sboxes[k])
        subsection(k)
        for table in ["LAT", "BCT", "DDT"]:
            print("{} {:10.5f} {:10.5f}".format(
                table,
                table_anomaly(s, table),
                table_negative_anomaly(s, table)
            ))
```

All the properties of the AES are extremely good, perhaps even optimal. Thus, the positive anomalies are essentially infinite (there are $256! \approx 2^{1684}$ 8-bit permutations, meaning an anomaly above 1684 doesn't really make sense). For the S-box of Fantomas, whose BCT properties are pretty bad, we observe the opposite. Finally, for the S-box of Kuznyechik, we can see that all the properties are too good to be random, especially the differential ones.


### Looking at the maximum value

`sboxU` provides tools to compute the expected value for the maximum coefficient in the LAT, DDT and BCT (based on [JMC:DaeRij07] and [AC:BonPerTia19]). Let's try them out in a new experiment.

```python
with Experiment("Testing distribution of maximum values"):
```

First, we set the global parameters. Statistical tests lose their relevance when the input size becomes too small so we fix `n=8`. Running a couple thousand tests is sufficient to gather enough information, hence `n_tested=2**13`.

```python
    n = 8
    n_tested = 2**13
```

#### LAT

```python
    section("LAT")
```

We first compute experimental data by generating `n_tested` random permutations. For each, we compute the linearity and increase a corresponding counter in a dictionary called `dis`. We then display the proportion of each linearity encountered.

```python
    subsection("Experimental")
    dis = defaultdict(int)
    for t in range(0, n_tested):
        dis[linearity(random_permutation_S_box(n))] += 1
    for k in sorted(dis.keys()):
        if dis[k] != 0:
            print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))
```

For comparison, we then generate the expected distribution using the relevant function (`expected_linearity_distribution_permutation`), and then display the result.

```python
    subsection("Theoretical")
    dis = expected_linearity_distribution_permutation(n, n)
    for k in sorted(dis.keys()):
        if dis[k] != 0.0:
            print("{:3d}: {:.4f}".format(k, dis[k]))

```


#### DDT

```python
    section("DDT")
```

Our experiment for the DDT is essentially the same as for the LAT and uses the same variables. We just need to replace `linearity` with `differential_uniformity` and `expected_linearity_distribution_permutation` with `expected_differential_uniformity_distribution_permutation`

```python
    subsection("Experimental")
    dis = defaultdict(int)
    for t in range(0, n_tested):
        dis[differential_uniformity(random_permutation_S_box(n))] += 1
    for k in sorted(dis.keys()):
        if dis[k] != 0:
            print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))
    
    subsection("Theoretical")
    dis = expected_differential_uniformity_distribution_permutation(n, n)
    for k in sorted(dis.keys()):
        if dis[k] != 0.0:
            print("{:3d}: {:.4f}".format(k, dis[k]))

```


#### BCT

```python
    section("BCT")
```

For the BCT, it is again the same except that we consider `boomerang_uniformity` with and `expected_boomerang_uniformity_distribution_permutations`.

```python
    subsection("Experimental")
    dis = defaultdict(int)
    for t in range(0, n_tested):
        dis[boomerang_uniformity(random_permutation_S_box(n))] += 1
    for k in sorted(dis.keys()):
        if dis[k] != 0:
            print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))
    
    subsection("Theoretical")
    dis = expected_boomerang_uniformity_distribution_permutation(n, n)
    for k in sorted(dis.keys()):
        if dis[k] != 0.0:
            print("{:3d}: {:.4f}".format(k, dis[k]))

```

#### Some Comments

As you can see by running this code, modeling the table coefficients like independent and identically distributed random variables works fairly well, though the results in the case of the expected linearity are not as good as the others. This is most likely due to the independance hypothesis being too heavy handed in this context.
