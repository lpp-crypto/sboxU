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
print(m_mat)
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

We can then easily check that all capture the same mathematical transformation.

```python
print("x\tm1(x)\tm2(x)\tmat*x")
for x in range(0, 2**3):
    print("{}\t{}\t{}\t{}".format(
    x, 
    m1(x),
    m2(x),
    from_bin(m_mat * vector(to_bin(x, 3)))))
```




## References


