# Basic Statistical Properties in the Binary Case



The corresponding source file is available online [on github](https://github.com/lpp-crypto/sboxU/blob/main/tests/statistics/tables.py).

## Preamble

Let's see how we can use `sboxU` to investigate the statistical properties of an S-box of $F_2^n$ in practice. To ease implementation, we will use the following packages.

```python
from collections import defaultdict
```

## Differential properties

The study of equations of the form $S(x+a)=S(x)+b$ is of crucial importance, for instance when investigating differential attacks [^C-BihSha91]. `sboxU` provides several utilities for this purpose.

First, let's pick a 6-bit permutation uniformly at random.

```python
s = random_permutation_S_box(6)
```

### Derivatives

First, it is possible to compute derivatives, i.e. given an S-box `s` to obtain the S-box corresponding to the vectorial boolean function $D_a s: x \mapsto s(x+a)+s(x)$, for any $a$. This is done using the `derivative` function. 

```python
D_1_s = derivative(s, 1)
pprint(D_1_s)
```

As a sanity check, we can verify that $D_a s(x) = D_a s(x + a)$, for all $x$.

```python
derivative_is_translation_invariant = True
for x in range(0, 2**s.get_output_length()):
    if D_1_s[x] != D_1_s[oplus(x, 1)]:
        fail("derivative should be identical on x and x+a for all x and a, but it isn't the case for x={}, a={}".format(x, a))
        derivative_is_translation_invariant = False
if derivative_is_translation_invariant:
    success("sanity check passed: the derivative on 1 is invariant under translation by 1")
```


### DDT

In general, it is convenient to compute the Difference Distribution Table (DDT). It is a table of integers of dimension $2^n \times 2^m$ such the entry `DDT[a][b]` is the number of solutions of the equation $s(x+a)+s(x)=b$. It is computed using the `ddt` function from `sboxU`.

```python
d = ddt(s)
```

Then, we can easily check the definition, reusing the derivatives `D_1_s` we computed above.

```python
ddt_row = [0 for x in s.input_space()]
for x in s.input_space():
    ddt_row[D_1_s[x]] += 1
if ddt_row == d[1]:
    success("The DDT row corresponding to input difference 1 is correct")
else:
    fail("Problem with the DDT")
```



## Linear properties


!TODO! talk about linear properties

### Boomerang properties

!TODO! talk about boomerang properties




## References

[^C-BihSha91]: Eli Biham and Adi Shamir. Differential cryptanalysis of des-like cryptosystems. In Alfred Menezes and Scott A. Vanstone, editors, Advances in Cryptology - CRYPTO '90, 10th Annual International Cryptology Conference, Santa Barbara, California, USA, August 11-15, 1990, Proceedings, Lecture Notes in Computer Science, 2–21. Springer, 1990. URL: https://doi.org/10.1007/3-540-38424-3_1, doi:10.1007/3-540-38424-3_1.

