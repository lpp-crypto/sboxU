# On Anomalies


The corresponding source file is available online [on github](https://github.com/lpp-crypto/sboxU/blob/main/tests/statistics/anomalies.py).

## Preamble

An "anomaly" is a concept that was introduced in [^AC-BonPerTia19] (though its principle was already present in [^C-BirPer15]). It can be positive or negative, and it quantifies how good the properties of a given S-box are compared to a permutation picked uniformly at random from the relevant set. If a an S-box has a high positive anomaly for a given property, then it means that its property is much better than should be expected from a random S-box. If instead its negative anomaly is high, then its property is much worse than should be expected.


Formally, consider a property $P$ defined for any S-box $S$, such that its outputs can be strictly ordered (for instance, if it's a number then we can simply compare two outputs).


Let's see how this notion behaves. We will need both the S-boxes from the literature that are bundled with SAGE, and the very useful `defaultdict` class.

```python
from sage.crypto.sboxes import sboxes
from collections import defaultdict
```


## Computing some anomalies

We will look at some S-boxes from the literature. To this end, we grab the `sboxes` dictionary shipped with SAGE which contains many of them.


### Skipjack

Let's start with the S-box of Skipjack and let's look at the anomalies corresponding to the main tables:

```python
skipjack = get_sbox(sboxes["Skipjack"])
for table in ["LAT", "BCT", "DDT"]:
    print("{} {:10.5f} {:10.5f}".format(
        table,
        table_anomaly(skipjack, table),
        table_negative_anomaly(skipjack, table)
    ))
```

As we can see, the positive anomaly for the LAT is very high (remember that it is a logarithm): the S-box has been somehow optimized to have good linear property (see [^C-BirPer15] and [^PhD-Perrin17] for more details). Another way this property can be seen is through the display of a graph of its absolute Walsh spectrum along with what is expected for a random permutation. This is only done using the following function, which will open a new window displaying this graph. The blue area also shows the discrepancy between the expected behaviour of a random S-box and that of Skipjack.

(you must close the graph window for the program to continue)

```python
interactive_distribution_comparison_lat(skipjack)
```

### Other S-boxes

We can look at these values for other S-boxes: those of the AES, Fantomas, and Kuznyechik.

```python
for k in ["AES", "Fantomas", "Kuznyechik"]:
    s = get_sbox(sboxes[k])
    print(k)
    for table in ["LAT", "BCT", "DDT"]:
        print("{} {:10.5f} {:10.5f}".format(
            table,
            table_anomaly(s, table),
            table_negative_anomaly(s, table)
        ))
```

All the properties of the AES are extremely good, perhaps even optimal. Thus, the positive anomalies are essentially infinite (there are $256! \approx 2^{1684}$ 8-bit permutations, meaning an anomaly above 1684 doesn't really make sense). For the S-box of Fantomas, whose BCT properties are pretty bad, we observe the opposite. Finally, for the S-box of Kuznyechik, we can see that all the properties are too good to be random, especially the differential ones.


## Looking at the maximum value

`sboxU` provides tools to compute the expected value for the maximum coefficient in the LAT, DDT and BCT (based on [^JMC-DaeRij07] and [^AC-BonPerTia19]). Let's try them out in a new experiment.

First, we set the global parameters. Statistical tests lose their relevance when the input size becomes too small so we fix `n=8`. Running a couple thousand tests is sufficient to gather enough information, hence `n_tested=2**13`.

```python
n = 8
n_tested = 2**11
```

### LAT 

#### Experimental 
We first compute experimental data by generating `n_tested` random permutations. For each, we compute the linearity and increase a corresponding counter in a dictionary called `dis`. We then display the proportion of each linearity encountered.

```python
dis = defaultdict(int)
for t in range(0, n_tested):
    dis[linearity(random_permutation_S_box(n))] += 1
for k in sorted(dis.keys()):
    if dis[k] != 0:
        print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))
```

#### Theoretical
For comparison, we then generate the expected distribution using the relevant function (`expected_linearity_distribution_permutation`), and then display the result.

```python
dis = expected_linearity_distribution_permutation(n, n)
for k in sorted(dis.keys()):
    if dis[k] != 0.0:
        print("{:3d}: {:.4f}".format(k, dis[k]))

```


### DDT

Our experiment for the DDT is essentially the same as for the LAT and uses the same variables. We just need to replace `linearity` with `differential_uniformity` and `expected_linearity_distribution_permutation` with `expected_differential_uniformity_distribution_permutation`

```python
print("Experimental")
dis = defaultdict(int)
for t in range(0, n_tested):
    dis[differential_uniformity(random_permutation_S_box(n))] += 1
for k in sorted(dis.keys()):
    if dis[k] != 0:
        print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))

print("Theoretical")
dis = expected_differential_uniformity_distribution_permutation(n, n)
for k in sorted(dis.keys()):
    if dis[k] != 0.0:
        print("{:3d}: {:.4f}".format(k, dis[k]))
```


### BCT

For the BCT, it is again the same except that we consider `boomerang_uniformity` with and `expected_boomerang_uniformity_distribution_permutations`.

```python
print("Experimental")
dis = defaultdict(int)
for t in range(0, n_tested):
    dis[boomerang_uniformity(random_permutation_S_box(n))] += 1
for k in sorted(dis.keys()):
    if dis[k] != 0:
        print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))

print("Theoretical")
dis = expected_boomerang_uniformity_distribution_permutation(n, n)
for k in sorted(dis.keys()):
    if dis[k] != 0.0:
        print("{:3d}: {:.4f}".format(k, dis[k]))
```

### Comments

As you can see by running this code, modeling the table coefficients like independent and identically distributed random variables works fairly well, though the results in the case of the expected linearity are not as good as the others. This is most likely due to the independance hypothesis being too heavy handed in this context.







## References

[^AC-BonPerTia19]: Xavier Bonnetain, Léo Perrin, and Shizhu Tian. Anomalies and vector space search: tools for s-box analysis. In Steven D. Galbraith and Shiho Moriai, editors, Advances in Cryptology - ASIACRYPT 2019 - 25th International Conference on the Theory and Application of Cryptology and Information Security, Kobe, Japan, December 8-12, 2019, Proceedings, Part I, volume 11921 of Lecture Notes in Computer Science, 196–223. Springer, 2019. URL: https://doi.org/10.1007/978-3-030-34578-5_8, doi:10.1007/978-3-030-34578-5_8.

[^C-BirPer15]: Alex Biryukov and Léo Perrin. On reverse-engineering s-boxes with hidden design criteria or structure. In Rosario Gennaro and Matthew Robshaw, editors, Advances in Cryptology - CRYPTO 2015 - 35th Annual Cryptology Conference, Santa Barbara, CA, USA, August 16-20, 2015, Proceedings, Part I, volume 9215 of Lecture Notes in Computer Science, 116–140. Springer, 2015. URL: https://doi.org/10.1007/978-3-662-47989-6_6, doi:10.1007/978-3-662-47989-6_6.

[^PhD-Perrin17]: Léo Perrin. Cryptanalysis, Reverse-Engineering and Design of Symmetric Cryptographic Algorithms. PhD thesis, University of Luxembourg, 2017. URL: http://orbilu.uni.lu/handle/10993/31195.

[^JMC-DaeRij07]: Joan Daemen and Vincent Rijmen. Probability distributions of correlation and differentials in block ciphers. J. Math. Cryptol., 1(3):221–242, 2007. URL: https://doi.org/10.1515/JMC.2007.011, doi:10.1515/JMC.2007.011.

