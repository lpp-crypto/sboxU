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


An `S_box` instance can be easily created using the LUT of the function investigated. For example, we generated the LUT of a random permutation of $F_{64}$: here is how to initialize an `S_box` instance using this LUT.

```python
s = get_sbox([57,29,43,15,16,12,10,59,37,63,32,23,9,24,31,1,33,
    26,39,38,18,7,58,49,27,54,62,52,4,53,5,25,47,
    11,40,0,50,48,8,42,35,13,3,36,2,44,56,6,28,
    34,17,14,45,30,61,21,51,19,46,60,41,20,55,22])
```

We can simply print `s` to see again the LUT, but it is more interesting to pretty print it using the `pprint` function. 

```python
print("basic list representation")
print(s)
print("pretty printing the same information")
pprint(s)
```

### Linear properties


### Differential properties


### Boomerang properties


### Polynomial representation




## Building an S-box

### The case of S-boxes from the literature

SUppose that you want to study the S-box of the AES (it has its own [wikipedia page](https://en.wikipedia.org/wiki/Rijndael_S-box). `sboxU` knows about the literature, and it is thus possible to write

```python
s = get_sbox("AES")
```

This will retrieve the LUT of the S-box of the AES from its internal database, and generate an object of the class `S_box` (which is assigned to `s`). 

```python
for test in [("AES", 4),
             ("Kuznyechik", 8)]:
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
