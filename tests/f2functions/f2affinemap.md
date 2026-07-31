# Affine Functions of F_2^n


The corresponding source file is available online [on github](https://github.com/lpp-crypto/sboxU/blob/main/tests/f2functions/f2affinemap.py).

## Preamble

Let $q=2^u$ be a power of two, and $F_q$ be the finite field with $q$ elements. `sboxU` provides fast tools to deal with affine functions mapping $F_q^n$ to $F_q^m$ (for now only in this case, meaning when the characteristic is 2), in particular the `F2AffineMap` class. They are implemented in C++, and the performance gain over plain SAGE matrices can be enormous. Since they are implemented in C++, you can also use them directly in C++ programs---but it is a topic for another time.

Under the hood, the `F2AffineMap` class defines matrices over $F_2$ only and stores both the offset and the image vectors using unsigned integers.

Let's see how it can be used in practice.


## Basic functionalities of the F2AffineMap class


### Construction 

SAGE has built-in support for matrices. For instance, the following builds the matrix corresponding to a linear bijection on 3 bits that leaves the bit of lowest weight unchanged, and swaps the two bits of highest weight.

```python
m_mat = Matrix(GF(2), 3, 3, [
    [1, 0, 0],
    [0, 0, 1],
    [0, 1, 0],
])
```

Such objects can be multiplied with each other, and multiplied with objects of the `vector` class.
The corresponding `F2AffineMap` can be built directly from this matrix.

```python
m1 = get_F2AffineMap(m_mat)
```

We can also build the same function by supplying the image vectors of the canonical basis.

```python
m2 = get_F2AffineMap([1, 4, 2])
```

As yet another method, we can supply the full lookup table of the functions. *Carefull* though: in this case, you need the input to be of type `S_box`, otherwise the constructor will assume that you are supplying images of the canonical basis! Even in this case, only the basis vectors are stored in the `F2AffineMap`.

```python
m3 = get_F2AffineMap(get_sbox([0, 1, 4, 5, 2, 3, 6, 7]))
```

We can then easily check that all capture the same mathematical transformation. Since the SAGE matrix operates on vectors of bits, we need to do some plumbing to work with it. Fortunately, the `sboxU` functions `from_bin` and `to_bin` functions simplify our lives. 

```python
print("x\tmat*x\tm1(x)\tm2(x)\tm3(x)")
for x in range(0, 2**3):
    y = from_bin(m_mat * vector(to_bin(x, 3)))
    row = "{}\t{}\t".format(x, y)
    for transformation in [m1, m2, m3]:
        y_prime = transformation(x)
        row += "{}\t".format(y_prime)
        if y_prime != y:
            fail("an F2AffineMap doesn't match the SAGE matrix")
    print(row)
```


### A bigger test: rank distribution

!TODO! write down the logic behind this part

```python
n_max = 10
parameters = [(7, 10), (10, 7), (10, 10)]
n_tested = 2**14
```

```python
def proba_full_rank(n, m):
    if m > n:
        m, n = n, m
    return float(prod(1 - 2**(-k) for k in range(n - m + 1, n + 1)))
```

```python
for n, m in parameters:
    counters = [0 for r in range(0, n_max+1)]
    for t in range(0, n_tested):
        L = rand_linear_function(n, m)
        counters[L.rank()] += 1
    row = "({:2d}, {:2d})".format(n, m)
    for c in counters:
        row += "\t{:5.3f}".format(float(c) / n_tested)
    print(row)
```

```python
    observed = float(counters[min(n, m)]) /  n_tested
    expected = proba_full_rank(n, m)
    diff = abs(expected - observed)
    if diff > 0.01:
        fail("mismatch between theory and practice: {}".format(diff))
    
```



## References


