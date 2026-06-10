# Dealing with S-boxes in F_2


The corresponding source file is available online [on github](https://github.com/lpp-crypto/sboxU/blob/main/tests/sbox/basicF2SBox.py).

## Preamble

Let $F_p$ be the finite field with $p$ elements, $p$ being a prime number (often, $p=2$, but not always!). `sboxU` is a library built to study *S-boxes*: we call *S-box* a function S satisfying the two following conditions:
- it maps $F_p^m$ to $F_p^n$, where both $n$ and $m$ are integers 
- its input size is small enough that that, in practice, it possible to precompute and store in memory all of its outputs.

There is a simple mapping $\lambda$ between the integers between 0 and $p^n-1$ and the elements of $F_p^n$. We then call the *lookup table* of S the sequence $S(\lambda(0)), ..., S(\lambda(p^n-1))$. 

In `sboxU`, the goal is to generate the LUT of a function in order to study it. This LUT is then wrapped in dedicated classes: `S_box` for S-boxes defined over vector spaces of $F_2$ (or, equivalently, over bitstrings), and `Fp_S_box` for other primes.

Let's see how this can be done in practice. 


## Basic functionalities of the S_box class


An `S_box` instance can be easily created using the LUT of the function investigated. For example,  here is how to initialize an `S_box` instance operating on $F_2^6$ that is simply the multiplication by 5 over the integers.

```python
N = 6
s = get_sbox([5*x % 2**N for x in range(0, 2**N)])
```

We can simply print `s` to see again the LUT, but it is more interesting to pretty print it using the `pprint` function. 

```python
print("basic list representation:")
print(s)
print("\npretty printing the same information:")
pprint(s)
```


### Queries

It is possible to query the content of the lookup of an S-box using the operator $[]$. The operator $()$ can also work but it offers additional functionalities (see the section on "Linear Casts").

For the all 0 vector, we have:
```python
if s[0] == 0:
    success("the output of s on the all zero bit string is indeed {}".format(s[0]))
else:
    fail("s[0] should be 0, something went wrong")
```

and for the vector `(1, 0, ..., 0)`, we need to use the integer 1 as the query. We get:

```python
if s[1] == 5:
    success("the output of s on (1,0,...,0) is indeed {}".format(s[0]))
else:
    fail("s[1] should be 5, something went wrong")
```

### Basic quantities


- input length
- output length
- input space
- input space size
- is_permutation

### Coordinates and Components

### Polynomial representations

- algegraic normal form
- univariate



## Building an S-box

### The case of S-boxes from the literature

SUppose that you want to study the S-box of the AES (it has its own [wikipedia page](https://en.wikipedia.org/wiki/Rijndael_S-box). `sboxU` knows about the literature, and it is thus possible to write

```python
s = get_sbox("AES")
```

This will retrieve the LUT of the S-box of the AES from its internal database, and generate an object of the class `S_box` (which is assigned to `s`). Thus, the following snippet grabs several S-boxes from the literature (namely, from the AES[^EC-DaeRij02], Ascon[^JC-DEMS21], and PRESENT[^CHES-BKLP+07].

```python
for test in [("AES", 4),
             ("PRESENT", 4),
             ("Ascon", 8)]:
    name, expected_uniformity = test
    u = differential_uniformity(get_sbox(name))
    if u == expected_uniformity:
        success("{} S-box differential uniformity is indeed {}".format(
            name,
            expected_uniformity
        ))
    else:
        fail("{} S-box differential uniformity should be {}, not {}".format(
            name,
            expected_uniformity,
            u
        ))

```


## Univariate polynomials


```python
field = GF(2**5)
g = field.gen()
X = PolynomialRing(field, "x").gen()
cube = get_sbox(X**3)
cube_plus = get_sbox(X**3 + g*X)
```


## Operations on S-boxes

### Addition

It is possible to add to S-boxes.

```python
diff = cube + cube_plus 
is_linear = (differential_uniformity(diff) == 2**5)
if is_linear:
    success("the sum of X^3 and X^3+gX is indeed an affine function")
else:
    fail("something went wrong, gX should be affine")
```


### Composition









## References

[^EC-DaeRij02]: Joan Daemen and Vincent Rijmen. AES and the wide trail design strategy. In Lars R. Knudsen, editor, Advances in Cryptology - EUROCRYPT 2002, International Conference on the Theory and Applications of Cryptographic Techniques, Amsterdam, The Netherlands, April 28 - May 2, 2002, Proceedings, Lecture Notes in Computer Science, 108–109. Springer, 2002. URL: https://doi.org/10.1007/3-540-46035-7_7, doi:10.1007/3-540-46035-7_7.

[^JC-DEMS21]: Christoph Dobraunig, Maria Eichlseder, Florian Mendel, and Martin Schläffer. Ascon v1.2: lightweight authenticated encryption and hashing. J. Cryptol., 34(3):33, 2021. URL: https://doi.org/10.1007/s00145-021-09398-9, doi:10.1007/S00145-021-09398-9.

[^CHES-BKLP+07]: Andrey Bogdanov, Lars R. Knudsen, Gregor Leander, Christof Paar, Axel Poschmann, Matthew J. B. Robshaw, Yannick Seurin, and C. Vikkelsoe. PRESENT: an ultra-lightweight block cipher. In Pascal Paillier and Ingrid Verbauwhede, editors, Cryptographic Hardware and Embedded Systems - CHES 2007, 9th International Workshop, Vienna, Austria, September 10-13, 2007, Proceedings, Lecture Notes in Computer Science, 450–466. Springer, 2007. URL: https://doi.org/10.1007/978-3-540-74735-2_31, doi:10.1007/978-3-540-74735-2_31.

