# S-boxes over F_p

The corresponding source file is available online [on github](https://github.com/lpp-crypto/sboxU/blob/main/tests/sbox/basicFpSBox.py).

## Preamble

This test file exercises the complete API of `S_box_fp`, the class used to represent vectorial functions over a prime field $\mathbb{F}_p$ with $p > 2$, that is, functions from $\mathbb{F}_p^n$ to $\mathbb{F}_p^m$ for some $n, m \in \mathbb{N}$.

```python
from sage.all import GF, PolynomialRing
from sboxU.core.sbox.cython_functions import S_box_fp
```

We define a main reference S-box `u` mapping $\mathbb{F}_3^2$ (so $m=n=2$ in this case) to itself.
It is one of the 9! bijections on this space and will be used throughout as a mean of example. We
represent its LUT as a list of lists of elements of $\mathbb{F}_3$, in
increasing order of the integer (little-endian base-3) representation of the
input.

```python
p = 3
Fp = GF(p)
raw_lut = [[0,1],[1,0],[0,2],[0,0],[2,0],[2,2],[1,2],[2,1],[1,1]]
lut = [[Fp(x) for x in row] for row in raw_lut]
```

We also define two helper functions used throughout the tests to compare
FpWord outputs. When talking about FpWord, we are refering to the type `list[Fp]` (the list of Fp-coordinates of the word):

```python
def fp_eq(a, b):
    """Return True if two FpWord lists are equal entry by entry."""
    return list(a) == list(b)

def fp_lut_eq(lut_a, lut_b):
    """Return True if two full LUTs are equal."""
    return all(fp_eq(a, b) for a, b in zip(lut_a, lut_b))
```


## Construction

We recall that the sbox that we want to construct is assumed to be from $\mathbb{F}_p^n$ to $\mathbb{F}_p^m$

### From a lookup table

The primary construction: pass a list of lists of $\mathbb{F}_p$ elements to
`get_sbox`. Each inner list is one output of the function in the canonical
ordering of the input space.

The (outer) list must then be of size $p^n$, while each inner list is of size $m$.

```python
u = get_sbox(lut)
if isinstance(u, S_box_fp):
    success("get_sbox(lut) correctly returned an S_box_fp instance")
else:
    fail("expected S_box_fp, got {}".format(type(u)))
```

### From multivariate polynomials

A list of multivariate polynomial ring over $\mathbb{F}_p$ can also be passed
directly. The LUT is then built by evaluating every polynomial on every point
of the input space.

Each polynomial must be lying in a polynomial field with $n$ variables, and the list must be made of $m$ variables.

Here we define a map $\mathbb{F}_3^2 \to \mathbb{F}_3^2$ whose first
coordinate is $x_1 x_2$ and whose second coordinate is $x_1 + x_2$.

```python
R = PolynomialRing(Fp, 2, "x")
x1, x2 = R.gens()
v = get_sbox([x1*x2, x1 + x2])
if isinstance(v, S_box_fp):
    success("get_sbox(polynomials) correctly returned an S_box_fp instance")
else:
    fail("expected S_box_fp from polynomial input, got {}".format(type(v)))
```

We can cross-check one entry: on input $(0,0)$ the expected output is
$(0 \cdot 0, 0+0) = (0,0)$.

```python
if fp_eq(v[[0,0]], [0,0]):
    success("v[(0,0)] = (0,0) as expected")
else:
    fail("v[(0,0)] = {} instead of (0,0)".format(list(v[[0,0]])))
```

And on input $(1,2)$ the expected output is $(1 \cdot 2, 1+2) = (2, 0)$.

```python
if fp_eq(v[[1,2]], [2,0]):
    success("v[(1,2)] = (2,0) as expected")
else:
    fail("v[(1,2)] = {} instead of (2,0)".format(list(v[[1,2]])))
```

### From a univariate polynomial

A univariate polynomial over $\mathbb{F}_p$ (or $\mathbb{F}_{p^n}$) can also
be passed to `get_sbox`. The LUT is built by evaluating the polynomial on
every element of the base field, using the polynomial basis to convert between
field elements and coordinate vectors.

For a polynomial over $\mathbb{F}_p$ (degree-1 case), the result has
`input_size = output_size = 1`.

```python
P3 = PolynomialRing(Fp, 'x')
xp = P3.gen()
sq_uni = get_sbox(xp**2)
if isinstance(sq_uni, S_box_fp):
    success("get_sbox(x**2 over GF(3)) returned an S_box_fp instance")
else:
    fail("expected S_box_fp from univariate polynomial, got {}".format(type(sq_uni)))
```

The squaring map over $\mathbb{F}_3$ is $0 \mapsto 0$, $1 \mapsto 1$, $2
\mapsto 1$, so the LUT should be `[[0],[1],[1]]`.

```python
expected_sq = [[0],[1],[1]]
got_sq = [list(y) for y in sq_uni.get_lut()]
if got_sq == expected_sq:
    success("squaring map x^2 over GF(3) has LUT [[0],[1],[1]]")
else:
    fail("squaring map LUT is {}, expected {}".format(got_sq, expected_sq))
```

We also check that the identity polynomial $x$ over $\mathbb{F}_3$ matches the
identity S-box, and that evaluating at a specific point is consistent with the
multivariate construction.

```python
id_uni = get_sbox(xp)
id_sp = [list(v) for v in id_uni.get_input_space()]
id_lut = [list(y) for y in id_uni.get_lut()]
if id_lut == id_sp:
    success("identity polynomial x over GF(3) gives the identity S-box")
else:
    fail("identity polynomial LUT is {}, expected {}".format(id_lut, id_sp))
```

For a polynomial over $\mathbb{F}_{p^n}$ the result has `input_size =
output_size = n`. Here we use $\mathbb{F}_9 = \mathbb{F}_{3^2}$ and the
identity map: every element should map to itself.

```python
K9 = GF(9, 'a')
P9 = PolynomialRing(K9, 'z')
z = P9.gen()
id9_uni = get_sbox(z)
sp9 = [list(v) for v in id9_uni.get_input_space()]
lut9 = [list(y) for y in id9_uni.get_lut()]
if lut9 == sp9 and id9_uni.get_input_size() == 2:
    success("identity polynomial z over GF(9) gives identity S-box on F_3^2")
else:
    fail("GF(9) identity polynomial: input_size={}, lut==input_space: {}".format(
        id9_uni.get_input_size(), lut9 == sp9))
```

### Bytes round-trip

The `to_bytes` method serializes any S-box to a byte string that `get_sbox`
can reconstruct. The leading byte is always `0x00` to distinguish this Fp
format from the F2 serialization.

This can be useful as a way to store on disk `S_box_fp` objects after a long computation. It is also used in the implementation of the `__hash()__` method for the class `S_box_fp`

```python
b = u.to_bytes()
if b[0] == 0:
    success("leading byte of to_bytes() is 0x00 as required by the Fp format")
else:
    fail("leading byte of to_bytes() is {}, expected 0".format(b[0]))
```

```python
u_reconstructed = get_sbox(b)
if isinstance(u_reconstructed, S_box_fp):
    success("get_sbox(to_bytes()) returns an S_box_fp")
else:
    fail("get_sbox(to_bytes()) returned {}, expected S_box_fp".format(
        type(u_reconstructed)))
```

```python
if u_reconstructed == u:
    success("round-trip to_bytes -> get_sbox is the identity")
else:
    fail("round-trip failed: reconstructed S-box differs from original")
```


## Basic properties

### Characteristic and sizes

```python
if u.get_p() == 3:
    success("get_p() returns 3")
else:
    fail("get_p() returned {}, expected 3".format(u.get_p()))
```

The input size is $n$.

```python
if u.get_input_size() == 2:
    success("get_input_size() returns 2")
else:
    fail("get_input_size() returned {}, expected 2".format(u.get_input_size()))
```
The output size is $m$.

```python
if u.get_output_size() == 2:
    success("get_output_size() returns 2")
else:
    fail("get_output_size() returned {}, expected 2".format(u.get_output_size()))
```
The size of the input size is $p^n$.

```python
if u.get_input_space_size() == 9:
    success("get_input_space_size() returns p^t = 9")
else:
    fail("get_input_space_size() returned {}, expected 9".format(u.get_input_space_size()))
```
The size of the output size is $p^m$

```python
if u.get_output_space_size() == 9:
    success("get_output_space_size() returns p^u = 9")
else:
    fail("get_output_space_size() returned {}, expected 9".format(u.get_output_space_size()))
```
The `len()` method is overriden, and returns the number of entries in the lookup table of the S_box.
 

```python
if len(u) == 9:
    success("len(u) returns 9")
else:
    fail("len(u) returned {}, expected 9".format(len(u)))
```

### Input and output spaces

The input space lists all $p^t = 9$ elements of $\mathbb{F}_3^2$ in
increasing order of their integer representation (little-endian base-3).

```python
expected_input_space = [
    [0,0],[1,0],[2,0],
    [0,1],[1,1],[2,1],
    [0,2],[1,2],[2,2]
]
actual = [list(x) for x in u.get_input_space()]
if actual == expected_input_space:
    success("get_input_space() enumerates all elements in correct (little-endian) order")
else:
    fail("get_input_space() returned {}, expected {}".format(actual, expected_input_space))
```

The output space has the same shape since input and output both live in
$\mathbb{F}_3^2$.

```python
if [list(x) for x in u.get_output_space()] == expected_input_space:
    success("get_output_space() matches the input space (same dimensions)")
else:
    fail("get_output_space() has unexpected content")
```

### Lookup table access

```python
lut_u = [list(x) for x in u.get_lut()]
if lut_u == raw_lut:
    success("get_lut() returns the original lookup table")
else:
    fail("get_lut() returned {}, expected {}".format(lut_u, raw_lut))
```

### String representation

`str(u)` should produce a list-of-lists representation of the lookup table in
the canonical format `[[a,b],[c,d],...]`.

```python
expected_str = "[[0,1],[1,0],[0,2],[0,0],[2,0],[2,2],[1,2],[2,1],[1,1]]"
if str(u) == expected_str:
    success("str(u) matches the expected representation")
else:
    fail("str(u) = '{}', expected '{}'".format(str(u), expected_str))
```


## Evaluation

### Direct lookup with __getitem__

`u[x]` takes an FpWord (a list or vector of integers/Fp elements) and returns the output.

```python
checks = [
    ([0,0], [0,1]),
    ([1,0], [1,0]),
    ([2,0], [0,2]),
    ([0,1], [0,0]),
    ([1,1], [2,0]),
    ([2,1], [2,2]),
    ([0,2], [1,2]),
    ([1,2], [2,1]),
    ([2,2], [1,1]),
]
all_ok = True
for x, expected in checks:
    if not fp_eq(u[x], expected):
        fail("u[{}] = {}, expected {}".format(x, list(u[x]), expected))
        all_ok = False
if all_ok:
    success("u[x] is correct for all 9 inputs")
```

### Iteration

Iterating over an S-box yields all output FpWords in input-space order.

```python
iterated = [list(y) for y in u]
if iterated == raw_lut:
    success("iterating over u yields all outputs in input-space order")
else:
    fail("iteration gave {}, expected {}".format(iterated, raw_lut))
```


## Equality, inequality and hash

Two S-boxes built from the same LUT must be equal.

```python
u2 = get_sbox(lut)
if u == u2:
    success("two S-boxes built from the same LUT are equal")
else:
    fail("u == u2 returned False, expected True")
```

An S-box is not equal to a different one.

```python
if u != v:
    success("u != v when the two S-boxes differ")
else:
    fail("u != v returned False for distinct S-boxes")
```

Hash must be consistent with equality.

```python
if hash(u) == hash(u2):
    success("equal S-boxes have the same hash")
else:
    fail("equal S-boxes have different hashes")
```


## Arithmetic operators

### Addition

Pointwise addition modulo $p$: `(u + v)[x] = u[x] + v[x]` coordinate-wise
mod 3. In particular, `u + u` doubles every output.

```python
uu = u + u
expected_double = [[(2*a) % 3 for a in row] for row in raw_lut]
got_double = [list(y) for y in uu]
if got_double == expected_double:
    success("u + u doubles every output coordinate mod 3")
else:
    fail("u + u gave {}, expected {}".format(got_double, expected_double))
```

### Subtraction

Pointwise subtraction modulo $p$: `(u - u)[x] = 0` for all $x$.

```python
u_minus_u = u - u
zero_lut = [[0, 0]] * 9
got_zero = [list(y) for y in u_minus_u]
if got_zero == zero_lut:
    success("u - u is the zero function")
else:
    fail("u - u gave {}, expected all zeros".format(got_zero))
```

Subtraction inverts addition: `(u + v) - v == u`.

```python
upv_minus_v = (u + v) - v
if upv_minus_v == u:
    success("(u + v) - v == u")
else:
    fail("(u + v) - v != u: got {}".format([list(y) for y in upv_minus_v]))
```

### Composition

Composition `u * w` applies `w` first, then `u`. Here we compose `u` with
itself: `(u * u)[x] = u(u(x))`.

```python
uu_comp = u * u
expected_comp = [list(u[list(row)]) for row in raw_lut]
got_comp = [list(y) for y in uu_comp]
if got_comp == expected_comp:
    success("(u * u)[x] == u(u(x)) for all x")
else:
    fail("composition gave {}, expected {}".format(got_comp, expected_comp))
```

### Power operator

`u ** 0` must be the identity on the input space.

```python
id_u = u ** 0
identity_lut = [list(x) for x in u.get_input_space()]
got_id = [list(y) for y in id_u]
if got_id == identity_lut:
    success("u ** 0 is the identity")
else:
    fail("u ** 0 gave {}, expected identity {}".format(got_id, identity_lut))
```

`u ** 2` must equal `u * u`.

```python
if (u ** 2) == (u * u):
    success("u ** 2 == u * u")
else:
    fail("u ** 2 != u * u")
```

`u ** -1` must equal `u.inverse()`.

```python
if (u ** -1) == u.inverse():
    success("u ** -1 == u.inverse()")
else:
    fail("u ** -1 != u.inverse()")
```


## Inversion

### Identity S-box

`S_box_fp.identity_S_box(t, p)` should return the identity map on
$\mathbb{F}_p^t$.

```python
id_sb = S_box_fp.identity_S_box(2, 3)
if fp_lut_eq(id_sb.get_lut(), u.get_input_space()):
    success("identity_S_box(2,3) maps every input to itself")
else:
    fail("identity_S_box(2,3) has wrong LUT: {}".format(
        [list(x) for x in id_sb.get_lut()]))
```

### Invertibility test

Our reference S-box `u` is a permutation of $\mathbb{F}_3^2$, so it must be
invertible. The polynomial-based `v` maps $(0,0) \to (0,0)$ and $(1,0) \to
(0,1)$, so we need to check whether it is also a bijection.

```python
if u.is_invertible():
    success("u.is_invertible() returns True for the bijection u")
else:
    fail("u.is_invertible() returned False for a known bijection")
```

### Inverse function

For a bijection $S$, we must have $S \circ S^{-1} = \mathrm{id}$ and
$S^{-1} \circ S = \mathrm{id}$.

```python
u_inv = u.inverse()
id_sb = S_box_fp.identity_S_box(2, 3)
```

```python
if (u * u_inv) == id_sb:
    success("u * u.inverse() is the identity")
else:
    fail("u * u.inverse() is not the identity")
```

```python
if (u_inv * u) == id_sb:
    success("u.inverse() * u is the identity")
else:
    fail("u.inverse() * u is not the identity")
```

Check one explicit entry: we have that
$u^{-1}([0,0]) = [0,1]$ (since $u([0,1]) = [0,0]$).

```python
if fp_eq(u_inv[[0,0]], [0,1]):
    success("u_inv[(0,0)] = (0,1) as expected from the LUT")
else:
    fail("u_inv[(0,0)] = {}, expected (0,1)".format(list(u_inv[[0,0]])))
```

## Coordinate functions and components

### Coordinate functions

`u.coordinate(i)` extracts the $i$-th scalar output coordinate: for every
$x$, `coordinate(i)[x] = [u[x][i]]`.

```python
coord_0 = u.coordinate(0)
expected_coord_0 = [[row[0]] for row in raw_lut]
got_coord_0 = [list(y) for y in coord_0]
if got_coord_0 == expected_coord_0:
    success("coordinate(0) correctly extracts the first output coordinate")
else:
    fail("coordinate(0) gave {}, expected {}".format(got_coord_0, expected_coord_0))
```

```python
coord_1 = u.coordinate(1)
expected_coord_1 = [[row[1]] for row in raw_lut]
got_coord_1 = [list(y) for y in coord_1]
if got_coord_1 == expected_coord_1:
    success("coordinate(1) correctly extracts the second output coordinate")
else:
    fail("coordinate(1) gave {}, expected {}".format(got_coord_1, expected_coord_1))
```

### Component functions

`u.component(a)` computes the scalar function $x \mapsto a \cdot S(x) \bmod p$
where $\cdot$ is the standard inner product over $\mathbb{F}_p$.

With $a = (1, 0)$, the component is the first coordinate function.

```python
comp_10 = u.component([1,0])
if comp_10 == coord_0:
    success("component([1,0]) equals coordinate(0)")
else:
    fail("component([1,0]) differs from coordinate(0): {} vs {}".format(
        [list(y) for y in comp_10], [list(y) for y in coord_0]))
```

With $a = (0, 1)$, the component is the second coordinate function.

```python
comp_01 = u.component([0,1])
if comp_01 == coord_1:
    success("component([0,1]) equals coordinate(1)")
else:
    fail("component([0,1]) differs from coordinate(1)")
```

For $a = (1, 1)$: the expected output at each input is the sum of both
coordinates mod 3.

```python
expected_comp_11 = [[(row[0] + row[1]) % 3] for row in raw_lut]
got_comp_11 = [list(y) for y in u.component([1,1])]
if got_comp_11 == expected_comp_11:
    success("component([1,1]) computes the sum of coordinates mod 3")
else:
    fail("component([1,1]) gave {}, expected {}".format(got_comp_11, expected_comp_11))
```

For $a = (2, 1)$: the expected output is $2 s_0 + s_1 \bmod 3$.

```python
expected_comp_21 = [[(2*row[0] + row[1]) % 3] for row in raw_lut]
got_comp_21 = [list(y) for y in u.component([2,1])]
if got_comp_21 == expected_comp_21:
    success("component([2,1]) computes 2*s0 + s1 mod 3")
else:
    fail("component([2,1]) gave {}, expected {}".format(got_comp_21, expected_comp_21))
```

The zero component $a = (0, 0)$ must be identically zero.

```python
comp_00 = u.component([0,0])
if [list(y) for y in comp_00] == [[0]] * 9:
    success("component([0,0]) is the zero function")
else:
    fail("component([0,0]) is not the zero function: {}".format(
        [list(y) for y in comp_00]))
```


## Derivative

The derivative of $S$ in direction $\delta$ is
$D_\delta(S)(x) = S(x + \delta) - S(x)$, where both the addition of $\delta$
and the coordinate-wise subtraction are done modulo $p$.

We first verify the trivial case: the zero derivative (with $\delta = 0$)
must be identically zero.

```python
zero_delta = [0, 0]
d_zero = u.derivative(zero_delta)
if [list(y) for y in d_zero] == [[0,0]] * 9:
    success("derivative in direction 0 is the zero function")
else:
    fail("derivative in direction 0 gave {}, expected all zeros".format(
        [list(y) for y in d_zero]))
```

For a non-trivial direction we compute the expected derivative independently
in Python and compare. The direction is $\delta = (1, 0)$.

```python
delta = [1, 0]
input_sp = [list(x) for x in u.get_input_space()]
input_size = u.get_input_size()
n = u.get_input_space_size()
powers = u.get_powers_in()

expected_deriv = []
for i in range(n):
    x = input_sp[i]
    x_plus_delta = [Fp(x[j]) + Fp(delta[j]) for j in range(input_size)]
    idx_xd = S_box_fp.vec_to_int(x_plus_delta, powers)
    out_xd = list(u.get_lut()[idx_xd])
    out_x  = list(u.get_lut()[i])
    expected_deriv.append([Fp(out_xd[j]) - Fp(out_x[j]) for j in range(u.get_output_size())])

got_deriv = [list(y) for y in u.derivative(delta)]
if got_deriv == expected_deriv:
    success("derivative in direction (1,0) matches the expected values")
else:
    fail("derivative in direction (1,0) gave {}, expected {}".format(
        got_deriv, expected_deriv))
```


## Serialization

The byte string returned by `to_bytes` encodes $p$, the output size, the
table length and the table body. We verify the header structure explicitly.

```python
import struct
bs = u.to_bytes()
# byte 0: format marker
# bytes 1..8: p (8-byte LE)
# bytes 9..16: output_size (8-byte LE)
# bytes 17..24: n = p^t (8-byte LE)
# byte 25: bytes_per_value
p_encoded      = struct.unpack_from('<q', bs, 1)[0]
out_encoded    = struct.unpack_from('<q', bs, 9)[0]
n_encoded      = struct.unpack_from('<q', bs, 17)[0]
bpv            = bs[25]
```

```python
if p_encoded == 3:
    success("p is correctly encoded as 3 in the header")
else:
    fail("p encoded as {}, expected 3".format(p_encoded))
```

```python
if out_encoded == 2:
    success("output_size is correctly encoded as 2 in the header")
else:
    fail("output_size encoded as {}, expected 2".format(out_encoded))
```

```python
if n_encoded == 9:
    success("n = p^t = 9 is correctly encoded in the header")
else:
    fail("n encoded as {}, expected 9".format(n_encoded))
```

```python
if bpv == 1:
    success("bytes_per_value is 1 (p=3 fits in a single byte)")
else:
    fail("bytes_per_value is {}, expected 1".format(bpv))
```

Finally, the total byte length must be $1 + 3 \times 8 + 1 + n \times
\text{output\_size} \times \text{bpv} = 26 + 9 \times 2 \times 1 = 44$.

```python
expected_len = 1 + 3*8 + 1 + 9 * 2 * 1
if len(bs) == expected_len:
    success("byte length of to_bytes() is {} as expected".format(expected_len))
else:
    fail("byte length is {}, expected {}".format(len(bs), expected_len))
```


## S-boxes over other primes

All tests so far used $p = 3$. We now confirm that the implementation is
agnostic to the specific prime.

### p = 5, bijection

We build the multiplication-by-2 map $x \mapsto 2x \bmod 5$ on
$\mathbb{F}_5$. Since $\gcd(2, 5) = 1$, it is a bijection.

```python
p5 = 5
Fp5 = GF(p5)
raw_lut_f = [[0],[2],[4],[1],[3]]
lut_f = [[Fp5(x) for x in row] for row in raw_lut_f]
f = get_sbox(lut_f)
if isinstance(f, S_box_fp):
    success("p=5 S-box is correctly an S_box_fp instance")
else:
    fail("p=5 S-box: expected S_box_fp, got {}".format(type(f)))
```

```python
if f.get_p() == 5 and f.get_input_size() == 1 and f.get_output_size() == 1 and len(f) == 5:
    success("p=5 S-box has correct characteristic, dimensions and length")
else:
    fail("p=5 S-box: wrong properties p={}, t={}, u={}, len={}".format(
        f.get_p(), f.get_input_size(), f.get_output_size(), len(f)))
```

```python
if fp_eq(f[[0]], [0]) and fp_eq(f[[1]], [2]) and fp_eq(f[[3]], [1]) and fp_eq(f[[4]], [3]):
    success("f[x] = 2x mod 5 is evaluated correctly")
else:
    fail("f evaluation incorrect")
```

```python
if f.is_invertible():
    success("multiplication by 2 in F_5 is correctly identified as a bijection")
else:
    fail("multiplication by 2 in F_5 should be a bijection")
```

Adding `f` to itself should give $4x \bmod 5$.

```python
ff = f + f
expected_4x = [[(4*x) % 5] for x in range(5)]
if [list(y) for y in ff] == expected_4x:
    success("f + f = 4x mod 5 for all x in F_5")
else:
    fail("f + f gave {}, expected {}".format([list(y) for y in ff], expected_4x))
```

`f - f` must be the zero function.

```python
if [list(y) for y in (f - f)] == [[0]] * 5:
    success("f - f is the zero function over F_5")
else:
    fail("f - f is not zero over F_5")
```

Bytes round-trip with $p = 5$.

```python
f_rt = get_sbox(f.to_bytes())
if isinstance(f_rt, S_box_fp) and f_rt == f:
    success("p=5 S-box survives to_bytes -> get_sbox round-trip")
else:
    fail("p=5 round-trip failed")
```

### p = 5, non-bijection

The squaring map $x \mapsto x^2 \bmod 5$ maps both $1$ and $4$ to $1$, and
both $2$ and $3$ to $4$, so it is not injective.

```python
raw_lut_g = [[0],[1],[4],[4],[1]]
lut_g = [[Fp5(x) for x in row] for row in raw_lut_g]
g = get_sbox(lut_g)
if not g.is_invertible():
    success("x^2 mod 5 is correctly identified as non-invertible")
else:
    fail("x^2 mod 5 should not be invertible")
```

### p = 7, round-trip

A cyclic shift by 1 over $\mathbb{F}_7$ is also a bijection. We use it
primarily to verify the bytes encoding for a larger prime.

```python
p7 = 7
Fp7 = GF(p7)
raw_lut_h7 = [[(x+1) % p7] for x in range(p7)]
lut_h7 = [[Fp7(x) for x in row] for row in raw_lut_h7]
h7 = get_sbox(lut_h7)
h7_rt = get_sbox(h7.to_bytes())
if isinstance(h7_rt, S_box_fp) and h7_rt == h7:
    success("p=7 S-box survives to_bytes -> get_sbox round-trip")
else:
    fail("p=7 round-trip failed")
```

The header of the byte string must encode $p = 7$.

```python
bs7 = h7.to_bytes()
p7_encoded = struct.unpack_from('<q', bs7, 1)[0]
if p7_encoded == 7:
    success("p=7 is correctly encoded in the bytes header")
else:
    fail("p=7 header encodes {}, expected 7".format(p7_encoded))
```


## Non-square S-boxes

### Injective: F_3^1 to F_3^2

The map $\phi \colon \mathbb{F}_3 \to \mathbb{F}_3^2$ defined by
$\phi(x) = (x, 2x \bmod 3)$ is injective but not surjective.

```python
raw_lut_phi = [[0,0],[1,2],[2,1]]
lut_phi = [[Fp(x) for x in row] for row in raw_lut_phi]
phi = get_sbox(lut_phi)
if phi.get_input_size() == 1 and phi.get_output_size() == 2:
    success("phi: F_3^1 -> F_3^2 has correct dimensions (t=1, u=2)")
else:
    fail("phi dimensions wrong: t={}, u={}".format(phi.get_input_size(), phi.get_output_size()))
```

```python
if len(phi) == 3 and phi.get_input_space_size() == 3 and phi.get_output_space_size() == 9:
    success("phi has 3 inputs and 9 possible outputs")
else:
    fail("phi wrong space sizes: len={}, in={}, out={}".format(
        len(phi), phi.get_input_space_size(), phi.get_output_space_size()))
```

```python
if fp_eq(phi[[0]], [0,0]) and fp_eq(phi[[1]], [1,2]) and fp_eq(phi[[2]], [2,1]):
    success("phi[x] = (x, 2x mod 3) evaluated correctly on all 3 inputs")
else:
    fail("phi evaluation failed: phi(0)={}, phi(1)={}, phi(2)={}".format(
        list(phi[[0]]), list(phi[[1]]), list(phi[[2]])))
```

```python
if not phi.is_invertible():
    success("phi is correctly non-invertible (input_size != output_size)")
else:
    fail("phi should not be invertible since input_size != output_size")
```

Bytes round-trip for non-square S-box.

```python
phi_rt = get_sbox(phi.to_bytes())
if isinstance(phi_rt, S_box_fp) and phi_rt == phi:
    success("non-square S-box (t=1, u=2) survives bytes round-trip")
else:
    fail("non-square S-box round-trip failed")
```

### Surjective: F_3^2 to F_3^1

The projection $\psi \colon \mathbb{F}_3^2 \to \mathbb{F}_3$ defined by
$\psi(x_0, x_1) = x_0 + x_1 \bmod 3$ is surjective but not injective.

```python
raw_lut_psi = [[0],[1],[2],[1],[2],[0],[2],[0],[1]]
lut_psi = [[Fp(x) for x in row] for row in raw_lut_psi]
psi = get_sbox(lut_psi)
if psi.get_input_size() == 2 and psi.get_output_size() == 1:
    success("psi: F_3^2 -> F_3^1 has correct dimensions (t=2, u=1)")
else:
    fail("psi dimensions wrong: t={}, u={}".format(psi.get_input_size(), psi.get_output_size()))
```

```python
if not psi.is_invertible():
    success("psi is correctly non-invertible (input_size != output_size)")
else:
    fail("psi should not be invertible")
```

Explicit checks: $\psi(0,0)=0$, $\psi(1,0)=1$, $\psi(2,1)=0$ (since $2+1=3 \equiv 0$).

```python
if fp_eq(psi[[0,0]], [0]) and fp_eq(psi[[1,0]], [1]) and fp_eq(psi[[2,1]], [0]):
    success("psi(x0,x1) = x0+x1 mod 3 evaluated correctly")
else:
    fail("psi evaluation failed")
```

Bytes round-trip for psi.

```python
psi_rt = get_sbox(psi.to_bytes())
if isinstance(psi_rt, S_box_fp) and psi_rt == psi:
    success("non-square S-box (t=2, u=1) survives bytes round-trip")
else:
    fail("non-square S-box (t=2, u=1) round-trip failed")
```

### Composing across different shapes

Composing $\psi \circ \phi \colon \mathbb{F}_3 \to \mathbb{F}_3$ applies
$\phi$ first (embedding into $\mathbb{F}_3^2$) then $\psi$ (projecting back).
We have $\psi(\phi(x)) = x + 2x = 3x \equiv 0 \pmod{3}$, so the result is
the constant zero function.

```python
psi_phi = psi * phi
if [list(y) for y in psi_phi] == [[0],[0],[0]]:
    success("psi * phi is the zero function (x + 2x = 0 mod 3)")
else:
    fail("psi * phi gave {}, expected all zeros".format([list(y) for y in psi_phi]))
```

Composing $u \circ \phi \colon \mathbb{F}_3 \to \mathbb{F}_3^2$ (embedding
into $\mathbb{F}_3^2$ then applying the bijection $u$). The output size of
$\phi$ is 2, which matches the input size of $u$.

```python
u_phi = u * phi
if u_phi.get_input_size() == 1 and u_phi.get_output_size() == 2:
    success("u * phi: F_3^1 -> F_3^2 has correct dimensions")
else:
    fail("u * phi has wrong dimensions: t={}, u={}".format(
        u_phi.get_input_size(), u_phi.get_output_size()))
```

```python
expected_u_phi = [list(u[list(row)]) for row in raw_lut_phi]
if [list(y) for y in u_phi] == expected_u_phi:
    success("(u * phi)[x] == u(phi(x)) for all x in F_3")
else:
    fail("u * phi gave {}, expected {}".format(
        [list(y) for y in u_phi], expected_u_phi))
```


## Error handling

The API must raise exceptions when the inputs are geometrically incompatible.

### Incompatible sizes in arithmetic

Adding two S-boxes with different lengths (here $p=3$ vs $p=5$, giving lengths
9 vs 5) must raise.

```python
try:
    _ = u + f
    fail("adding S-boxes with different lengths should raise")
except Exception:
    success("u + f raises when lengths differ (9 vs 5)")
```

```python
try:
    _ = u - f
    fail("subtracting S-boxes with different lengths should raise")
except Exception:
    success("u - f raises when lengths differ")
```

### Incompatible sizes in composition

`u * psi`: the output size of $\psi$ is 1, but the input size of $u$ is 2.

```python
try:
    _ = u * psi
    fail("composing u (input_size=2) with psi (output_size=1) should raise")
except Exception:
    success("u * psi raises: output_size(psi)=1 != input_size(u)=2")
```

`phi * u`: the output size of $u$ is 2, but the input size of $\phi$ is 1.

```python
try:
    _ = phi * u
    fail("composing phi (input_size=1) with u (output_size=2) should raise")
except Exception:
    success("phi * u raises: output_size(u)=2 != input_size(phi)=1")
```

### Inverting a non-bijection

```python
try:
    _ = g.inverse()
    fail("g.inverse() should raise for the non-bijective x^2 mod 5")
except Exception:
    success("g.inverse() correctly raises for non-bijective g")
```

```python
try:
    _ = g ** -1
    fail("g ** -1 should raise for non-bijective g")
except Exception:
    success("g ** -1 correctly raises for non-bijective g")
```

### component with wrong vector length

`u` has output size 2; passing a vector of length 3 must raise.

```python
try:
    _ = u.component([1, 0, 0])
    fail("component([1,0,0]) on output_size=2 S-box should raise")
except Exception:
    success("component([1,0,0]) on output_size=2 S-box correctly raises")
```

### coordinate out of range

`u` has output size 2; valid indices are 0 and 1.

```python
try:
    _ = u.coordinate(2)
    fail("coordinate(2) on output_size=2 S-box should raise")
except Exception:
    success("coordinate(2) on output_size=2 S-box correctly raises")
```


## Mathematical properties

### Commutativity and associativity of addition

```python
if (u + v) == (v + u):
    success("addition is commutative: u + v == v + u")
else:
    fail("addition is not commutative!")
```

Associativity: $(u + v) + w = u + (v + w)$. We use $w = u$ as a third operand.

```python
if ((u + v) + u) == (u + (v + u)):
    success("addition is associative: (u + v) + u == u + (v + u)")
else:
    fail("addition is not associative!")
```

Zero is the additive identity: $u + 0 = u$ where $0$ is $u - u$.

```python
zero_sb = u - u
if (u + zero_sb) == u:
    success("u + zero == u: zero is the additive identity")
else:
    fail("u + zero != u")
```

### Identity law for composition

Composing any S-box with the identity on either side gives back the same S-box.

```python
if (u * id_sb) == u and (id_sb * u) == u:
    success("u * Id == u and Id * u == u")
else:
    fail("identity law for composition failed")
```

### Higher power

$u^3 = u \circ u \circ u$.

```python
if (u ** 3) == (u * u * u):
    success("u ** 3 == u * u * u")
else:
    fail("u ** 3 != u * u * u")
```

### Inverse undoes iteration

$(u^2)^{-1} = (u^{-1})^2$, since $(u^{-1})^2 \circ u^2 = \mathrm{id}$.

```python
if (u ** 2).inverse() == (u.inverse() ** 2):
    success("(u ** 2).inverse() == u.inverse() ** 2")
else:
    fail("(u ** 2).inverse() != u.inverse() ** 2")
```

### Derivative in multiple directions

We verify the derivative against our reference Python implementation for two
additional directions, $\delta = (0, 1)$ and $\delta = (2, 2)$.

```python
def ref_deriv(sbox, delta):
    Fp_s = GF(sbox.get_p())
    lut_s = [list(y) for y in sbox.get_lut()]
    sp = [list(x) for x in sbox.get_input_space()]
    powers = sbox.get_powers_in()
    in_dim = sbox.get_input_size()
    out_dim = sbox.get_output_size()
    result = []
    for i, x in enumerate(sp):
        xd = [Fp_s(x[j]) + Fp_s(delta[j]) for j in range(in_dim)]
        idx_xd = S_box_fp.vec_to_int(xd, powers)
        result.append([Fp_s(lut_s[idx_xd][j]) - Fp_s(lut_s[i][j]) for j in range(out_dim)])
    return result
```

```python
delta_01 = [0, 1]
got_01 = [list(y) for y in u.derivative(delta_01)]
exp_01 = ref_deriv(u, delta_01)
if got_01 == exp_01:
    success("derivative in direction (0,1) matches reference")
else:
    fail("derivative in direction (0,1): got {}, expected {}".format(got_01, exp_01))
```

```python
delta_22 = [2, 2]
got_22 = [list(y) for y in u.derivative(delta_22)]
exp_22 = ref_deriv(u, delta_22)
if got_22 == exp_22:
    success("derivative in direction (2,2) matches reference")
else:
    fail("derivative in direction (2,2): got {}, expected {}".format(got_22, exp_22))
```

### Derivative of the identity gives delta

For the identity map, $D_\delta(\mathrm{Id})(x) = \mathrm{Id}(x+\delta) -
\mathrm{Id}(x) = \delta$ for all $x$. So the derivative of the identity is
the constant S-box that always outputs $\delta$.

```python
delta_12 = [1, 2]
d_id = id_sb.derivative(delta_12)
if [list(y) for y in d_id] == [delta_12] * 9:
    success("derivative of identity is the constant delta=(1,2) for all x")
else:
    fail("derivative of identity gave {}, expected constant {}".format(
        [list(y) for y in d_id], [delta_12]*9))
```

### Linearity of the derivative in S

$D_\delta(S + T) = D_\delta(S) + D_\delta(T)$ for any two S-boxes $S$, $T$ of
matching dimensions and any direction $\delta$.

```python
delta_lin = [1, 0]
du = u.derivative(delta_lin)
dv = v.derivative(delta_lin)
dupv = (u + v).derivative(delta_lin)
if dupv == (du + dv):
    success("D_delta(u+v) == D_delta(u) + D_delta(v): derivative is linear in S")
else:
    fail("derivative linearity failed")
```

### Derivative of the zero function is zero

If $S$ is identically zero, $D_\delta(S)(x) = 0 - 0 = 0$ for all $x$ and
$\delta$.

```python
d_zero_fn = zero_sb.derivative([2, 1])
if [list(y) for y in d_zero_fn] == [[0,0]] * 9:
    success("derivative of the zero function is zero for any direction")
else:
    fail("derivative of zero function is not zero: {}".format(
        [list(y) for y in d_zero_fn]))
```

### Derivation is commutative

When deriving two times $D_a D_b(S)(x) = D_b D_a(S)(x)$ must hold for any two
directions $a$, $b$.

```python
da = [1, 0]
db = [0, 1]
d_a_then_b = u.derivative(da).derivative(db)
d_b_then_a = u.derivative(db).derivative(da)
if d_a_then_b == d_b_then_a:
    success("second derivative is symmetric: D_a D_b(S) == D_b D_a(S)")
else:
    fail("second derivative not symmetric!")
```

### Coordinate of derivative equals derivative of coordinate

$D_\delta(S)_i = D_\delta(S_i)$, where $S_i$ is the $i$-th coordinate of $S$.

```python
delta_cd = [1, 0]
deriv_coord_0 = u.derivative(delta_cd).coordinate(0)
coord_0_deriv = u.coordinate(0).derivative(delta_cd)
if deriv_coord_0 == coord_0_deriv:
    success("D_delta(u)_0 == D_delta(u_0): derivative commutes with coordinate extraction")
else:
    fail("D_delta(u)_0 != D_delta(u_0)")
```

### Equality is reflexive, symmetric, and consistent with iteration

```python
u3 = get_sbox(lut)
if u == u3 and u3 == u:
    success("equality is symmetric")
else:
    fail("equality is not symmetric")
```

```python
if u == u:
    success("equality is reflexive")
else:
    fail("equality is not reflexive!")
```

Different S-boxes must not be equal even if they share $p$ and dimensions.

```python
u_shifted = u + u - u  # should equal u
if u_shifted == u:
    success("u + u - u == u: arithmetic round-trips to equality")
else:
    fail("u + u - u != u")
```

```python
if (u + u) != u:
    success("u + u != u for non-trivial u")
else:
    fail("u + u == u, which should not happen for non-zero u")
```

## DDT properties

The *Difference Distribution Table* (DDT) of an $\mathbb{F}_p$ S-box $S : \mathbb{F}_p^n \to \mathbb{F}_p^m$ is a two-dimensional array
$$D[\delta][\gamma] = \#\{x \in \mathbb{F}_p^n : S(x+\delta) - S(x) = \gamma\}$$
where $+$ and $-$ are coordinate-wise modular arithmetic.

```python
from sboxU.statistics import fp_ddt, fp_differential_spectrum, fp_differential_uniformity
```

We use two simple reference S-boxes over $\mathbb{F}_3$:

```python
id3 = get_sbox([[Fp(0)],[Fp(1)],[Fp(2)]])
sq3 = get_sbox([[Fp(0)],[Fp(1)],[Fp(1)]])
```

`id3` is the identity on $\mathbb{F}_3$, and `sq3` is the squaring map $x \mapsto x^2 \bmod 3$, which takes values $[0,1,1]$ and is a *perfect nonlinear* (PN) function.

### Row sums equal input-space size

Every row of the DDT sums to the size of the input space $p^n$, since for every fixed $\delta$ the differences $S(x+\delta)-S(x)$ partition the $p^n$ inputs.

```python
ddt_id = fp_ddt(id3)
if all(sum(row) == 3 for row in ddt_id):
    success("every DDT row of id3 sums to 3")
else:
    fail("some DDT row of id3 has wrong sum: {}".format([sum(r) for r in ddt_id]))
```

### Delta=0 row

For $\delta = 0$ we have $S(x+0) - S(x) = 0$ for all $x$, so $D[0][0] = p^n$ and all other entries in row 0 are zero.

```python
if ddt_id[0] == [3, 0, 0]:
    success("DDT row 0 of id3 is [3,0,0]")
else:
    fail("DDT row 0 of id3 is {}, expected [3,0,0]".format(ddt_id[0]))
```

### Identity S-box rows

For the identity $S(x) = x$, the output difference equals the input difference: $D[\delta][\gamma] = p^n$ if $\gamma = \delta$ and $0$ otherwise.

```python
if ddt_id[1] == [0, 3, 0]:
    success("DDT row 1 of id3 is [0,3,0]")
else:
    fail("DDT row 1 of id3 is {}, expected [0,3,0]".format(ddt_id[1]))
```

```python
if ddt_id[2] == [0, 0, 3]:
    success("DDT row 2 of id3 is [0,0,3]")
else:
    fail("DDT row 2 of id3 is {}, expected [0,0,3]".format(ddt_id[2]))
```

### PN function: uniform DDT rows

A perfect nonlinear function has all DDT entries equal to $p^{n-m}$ for $\delta \ne 0$.
For `sq3` ($n=m=1$, $p=3$) this means every non-zero row should be $[1,1,1]$.

```python
ddt_sq = fp_ddt(sq3)
if ddt_sq[1] == [1, 1, 1] and ddt_sq[2] == [1, 1, 1]:
    success("DDT of squaring map has uniform rows [1,1,1] for nonzero delta")
else:
    fail("DDT of squaring map: row1={}, row2={}".format(ddt_sq[1], ddt_sq[2]))
```

### Differential spectrum

The differential spectrum counts, for each coefficient value $k$, the number of pairs $(\delta, \gamma)$ with $\delta \ne 0$ and $D[\delta][\gamma] = k$.

```python
sp_id = fp_differential_spectrum(id3)
from sboxU.core.spectrum.cython_functions import Spectrum
if isinstance(sp_id, Spectrum):
    success("fp_differential_spectrum returns a Spectrum instance")
else:
    fail("fp_differential_spectrum returned {}, expected Spectrum".format(type(sp_id)))
```

For `id3`, the non-zero delta rows each have exactly one entry equal to 3 (and two entries equal to 0).
So the spectrum has $6$ zeros and $2$ threes among the $2 \times 3 = 6$ (delta, gamma) pairs with delta != 0.

```python
if sp_id[3] == 2:
    success("fp_differential_spectrum(id3)[3] == 2")
else:
    fail("fp_differential_spectrum(id3)[3] == {}, expected 2".format(sp_id[3]))
```

```python
if sp_id[0] == 4:
    success("fp_differential_spectrum(id3)[0] == 4")
else:
    fail("fp_differential_spectrum(id3)[0] == {}, expected 4".format(sp_id[0]))
```

### Differential uniformity

The differential uniformity is the maximum DDT coefficient over all $\delta \ne 0$.

```python
if fp_differential_uniformity(id3) == 3:
    success("differential uniformity of id3 is 3")
else:
    fail("differential uniformity of id3 is {}, expected 3".format(fp_differential_uniformity(id3)))
```

```python
if fp_differential_uniformity(sq3) == 1:
    success("differential uniformity of squaring map (PN) is 1")
else:
    fail("differential uniformity of sq3 is {}, expected 1".format(fp_differential_uniformity(sq3)))
```
