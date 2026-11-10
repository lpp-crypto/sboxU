# Basic Operations on F_2 Vectors

The corresponding source file is available [online]()

## Preamble

`sboxU` represents elements of $(F_2)^n$ (i.e. bit vectors of length $n$) as plain integers, which the codebase calls `BinWord`. Bit $i$ of the integer `x` (i.e. `(x >> i) & 1`) corresponds to the $i$-th coordinate of the vector, with bit $0$ being the *least significant bit*. This is the convention used throughout `sboxU.core.f2functions`, in particular for the functions we look at here: `lsb`, `msb`, `to_bin`, `from_bin`, `scal_prod` and `hamming_weight`. Under the hood, all of them are thin Cython wrappers around small C++ intrinsics, so they are essentially free to call.

To make the demonstrations below easier to read, we define a small helper that prints the binary decomposition of an integer in the usual, human-readable, most-significant-bit-first order.

```python
def pretty_bin(x, n):
    return "".join(str(b) for b in reversed(to_bin(x, n)))
    
prng = InsecurePRNG(b"seed")
```

## Converting between integers and bit lists: `to_bin` and `from_bin`

`to_bin(x, n)` turns the `BinWord` `x` into a list of `n` bits, indexed from the least significant one (position `0`) to the most significant one (position `n-1`). `from_bin` performs the inverse operation, turning such a list back into a `BinWord`.

```python
n = 8
x = 0b01101001
print("x           =", pretty_bin(x, n))
print("to_bin(x,n) =", to_bin(x, n))
```

As expected, `to_bin(x, n)[0]` is the least significant bit of `x`, and `to_bin(x, n)[n-1]` is the most significant one (within the `n`-bit window we asked for).

The two functions are inverse of each other, which we can check on a handful of random values:

```python
from sage.all import randint

ok = True
for i in range(1000):
    y = prng(0, 2**n)
    if from_bin(to_bin(y, n)) != y:
        ok = False
if not ok:
    fail("to_bin/from_bin round-trip")
```

## Extremal bits: `lsb` and `msb`

`lsb(x)` and `msb(x)` return the position of the least and most significant bits of `x` that are set to `1`, i.e. the smallest and largest `i` such that `(x >> i) & 1 == 1`. By convention, both return `0` when `x` is `0`.

```python
for x in [0b00000001, 0b00001000, 0b01101001, 0]:
    print("x = {:>8s}   lsb(x) = {}   msb(x) = {}".format(
        pretty_bin(x, n), lsb(x), msb(x)
    ))
```

A useful identity, which we can check directly, is that `x >> msb(x)` is always exactly `1`, as long as `x` is not `0`:

```python
for _ in range(10):
    y = randint(1, 2**n - 1)
    assert (y >> msb(y)) == 1
print("x >> msb(x) == 1 for all tested non-zero x")
```

## Counting set bits: `hamming_weight`

`hamming_weight(x)` returns the number of bits equal to `1` in the binary decomposition of `x`, i.e. the Hamming weight of the corresponding vector in $(F_2)^n$.

```python
for x in [0, 1, 0b01101001, 2**n - 1]:
    print("x = {:>8s}   hamming_weight(x) = {}".format(pretty_bin(x, n), hamming_weight(x)))
```

Note that this matches simply summing up the entries returned by `to_bin`:

```python
for _ in range(10):
    y = randint(0, 2**n - 1)
    assert hamming_weight(y) == sum(to_bin(y, n))
print("hamming_weight(x) == sum(to_bin(x,n)) for all tested x")
```

## The canonical scalar product: `scal_prod`

`scal_prod(x, y)` computes the canonical scalar product of `x` and `y` over $F_2$, i.e. the parity (sum modulo 2) of the bitwise AND of `x` and `y`. This is one of the most commonly used building blocks in the rest of `sboxU`, in particular when computing components of S-boxes or Walsh transforms.

```python
x = 0b01101001
y = 0b01000011
print("x            =", pretty_bin(x, n))
print("y            =", pretty_bin(y, n))
print("x AND y      =", pretty_bin(x & y, n))
print("scal_prod(x,y) =", scal_prod(x, y))
```

As the name suggests, this is exactly the parity of the Hamming weight of `x & y`, which ties this section back to `hamming_weight`:

```python
for i in range(10):
    a = prng(0, 2**n)
    b = prng(0, 2**n)
    assert scal_prod(a, b) == hamming_weight(a & b) % 2
print("scal_prod(x,y) == hamming_weight(x & y) % 2 for all tested x,y")
```

## Comments

These six functions look elementary, but they are exactly the primitives that the rest of `sboxU` builds on: `to_bin`/`from_bin` give the bridge between `BinWord` integers and explicit vectors, `lsb`/`msb`/`hamming_weight` are the basic tools used when enumerating or characterizing vectors (e.g. cosets, weight classes), and `scal_prod` is the building block behind components, coordinates, and Walsh/Fourier-style transforms used throughout the rest of the library.
