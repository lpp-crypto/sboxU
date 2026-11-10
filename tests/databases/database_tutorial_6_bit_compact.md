# Tutorial - Generating compact 6-bit APN database

The goal of this tutorial is to show how to build a database with the class APNFunctions_compact. To do so, we are going to build the database of the known APN functions in dimension 6.

Building the database will be made in two steps. We first insert the representatives of the known CCZ-class: there are 13 functions from the Banff list, source 0 in the source dictionary, and 1 function that is non ccz-quadratic from Edel and Pott, source 1 in the source dictionary. In the second step, we add the ccz classes of these functions.

Executing this script will take few minutes due to the non ccz-quadratic class, as the other are treated in few seconds


## Preamble

```python

import os
from sboxU.core import get_sbox
from sage.all import GF, PolynomialRing

N_QUADRATIC = 13  # indices 0..12 in ccz_class_representatives

field = GF(2**6)
X = PolynomialRing(field, "X").gen()
g = field.gen()

ccz_class_representatives = [
    # Banff list for quadratic functions
    get_sbox(X**3),
    get_sbox(X**3 + g**11*X**6 + g*X**9),
    get_sbox(g*X**5 + X**9 + g**4*X**17 + g*X**18 + g**4*X**20 + g*X**24 + g**4*X**34 + g*X**40),
    get_sbox(g**7*X**3 + X**5 + g**3*X**9 + g**4*X**10 + X**17 + g**6*X**18),
    get_sbox(X**3 + g*X**24 + X**10), # 4 <-- KIM
    get_sbox(X**3 + g**17*(X**17 + X**18 + X**20 + X**24)),
    get_sbox(X**3 + g**11*X**5 + g**13*X**9 + X**17 + g**11*X**33 + X**48),
    get_sbox(g**25*X**5 + X**9 + g**38*X**12 + g**25*X**18 + g**25*X**36),
    get_sbox(g**40*X**5 + g**10*X**6 + g**62*X**20 + g**35*X**33 + g**15*X**34 + g**29*X**48),
    get_sbox(g**34*X**6 + g**52*X**9 + g**48*X**12 + g**6*X**20 + g**9*X**33 + g**23*X**34 + g**25*X**40),
    get_sbox(X**9 + g**4*(X**10 + X**18 ) + g**9*(X**12 + X**20 + X**40 )),
    get_sbox(g**52*X**3 + g**47*X**5 + g*X**6 + g**9*X**9 + g**44*X**12 + g**47*X**33 + g**10*X**34 + g**33*X**40),
    get_sbox(g*(X**6 + X**10 + X**24 + X**33) + X**9 + g**4*X**17),
    # cubic function from Edel and Pott
    get_sbox([0, 0, 0, 8, 0, 26, 40, 58, 0, 33, 10, 35, 12, 55, 46, 29, 0, 11, 12, 15, 4, 21, 32, 57, 20, 62, 18, 48, 28, 44, 50, 10, 0, 6, 18, 28, 10, 22, 48, 36, 8, 47, 16, 63, 14, 51, 62, 11, 5, 24, 27, 14, 11, 12, 61, 50, 25, 37, 13, 57, 27, 61, 39, 9])
]


```


## Setup

We simply remove the database if it does exist, and write over it if it exists.
```python

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "compact6.db")
if os.path.exists(DB_PATH):
    print("Removing existing database: {}".format(DB_PATH))
    os.remove(DB_PATH)
print("Database will be written to: {}".format(DB_PATH))
```


## Insert and populate database

We specify the path and  the dimension as a parameter
```python
with APNFunctions_compact(DB_PATH, n=6) as db:

```
We then  add the quadratic functions

```python
    section("Inserting quadratic CCZ-class representatives")
    print("(source 0 = Banff list)")
    quad_entry_ids = []
    for i in range(N_QUADRATIC):
        eid = db.insert_quadratic(ccz_class_representatives[i], source=0)
        quad_entry_ids.append(eid)
        print("  [{:2d}/{}] ccz_id {:2d} inserted → entry id {}".format(i + 1, N_QUADRATIC, i, eid))
    success("Inserted {} quadratic representatives".format(N_QUADRATIC))

```
And populate the classes
```python
    section("Populating quadratic CCZ classes")
    for i, eid in enumerate(quad_entry_ids):
        subsection("ccz_id {:2d}  (entry id {})".format(i, eid))
        inserted = db.populate_quadratic_ccz_class(eid)
        print("  {} EA representatives added".format(len(inserted)))
    success("All {} quadratic CCZ classes populated".format(N_QUADRATIC))

```
We now repeat the same process for the non-quadratic class, this step is also fast since it is a simple insertion.


```python
    section("Inserting non-quadratic CCZ-class representatives")
    print("(source 1 = Edel and Pott)")
    quad_entry_ids = []

    eid = db.insert_non_quadratic(ccz_class_representatives[13], source=1)
    quad_entry_ids.append(eid)
    print("  [1/1] non-quadratic representative inserted → entry id {}".format(eid))
    success("Inserted 1 non-quadratic representative")

```

We now populate the classes, which can take few minutes since function 13 is not ccz-quadratic.

```python
    section("Populating non-quadratic CCZ classes")
    print("This section should take few minutes")
    for i, eid in enumerate(quad_entry_ids):
        subsection("ccz_id {:2d}  (entry id {})".format(i, eid))
        inserted = db.populate_non_quadratic_ccz_class(eid)
        print("  {} EA representatives added".format(len(inserted)))
    success("All {} non-quadratic CCZ classes populated".format(1))

```