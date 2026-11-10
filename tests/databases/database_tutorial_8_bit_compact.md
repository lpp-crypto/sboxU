# Tutorial - Generating the first compact 8-bit APN database

The goal of this tutorial is to show how to build a database with the class APNFunctions_compact for 8-bit functions. To do so, we are going to batch-insert the known quadratic CCZ-class representatives from five different sources.

Unlike the 6-bit case, only insertion is performed here: CCZ classes are NOT populated with their further EA-class representatives, since doing so for the ~33k quadratic representatives inserted here would be prohibitively expensive.

Each source is tagged with its own bibliographic source index, taken from the biblio_8 dictionary in sboxU.apn.biblio_dict:
- source 5, Weng, Tan & Gong (2013): 10 functions
- source 1, Yu, Wang & Li (2014) - first QAM: 8179 functions
- source 2, Yu & Perrin (2022) - second QAM: 5412 functions
- source 3, Beierle & Leander (2022): 12923 functions
- source 4, Beierle, Leander & Perrin (2022): 6368 functions

Each group is inserted with a single call to batch_insert_quadratic_function, with a comment describing where the functions came from. 
We also update the database's journal as we insert functions to keep track of the modifications.


## Preamble

```python

import os
import sys
from sboxU.apn import APNFunctions_compact
# These functions are already stored in SboxU
from sboxU.scripts.apnDB.reprs8 import (
    all_WenTanGon, first_QAMs, second_QAMs, all_BeiLea, all_BLP22,
)

# (function group getter, human label, biblio_8 source index)
APNS = [
    (all_WenTanGon, "Weng, Tan & Gong (2013)",           5),
    (first_QAMs,    "Yu, Wang & Li (2014) - first QAM",  1),
    (second_QAMs,   "Yu & Perrin (2022) - second QAM",   2),
    (all_BeiLea,    "Beierle & Leander (2022)",           3),
    (all_BLP22,     "Beierle, Leander & Perrin (2022)",   4),
]

```


## Setup

We simply remove the database if it does exist, and write over it if it exists.
```python

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "compact8.db")
if os.path.exists(DB_PATH):
    print("Removing existing database: {}".format(DB_PATH))
    os.remove(DB_PATH)
print("Database will be written to: {}".format(DB_PATH))
```


## Insert database

We specify the path and the dimension as a parameter
```python
with APNFunctions_compact(DB_PATH, n=8) as db:

```
We then batch-insert each source group in turn, tagging every batch with a comment naming where it comes from

```python
    section("Batch-inserting quadratic CCZ-class representatives")
    for group_fn, label, source in APNS:
        subsection(label)
        functions = group_fn()
        print("  {} functions (source {})".format(len(functions), source))
        ids = db.batch_insert_quadratic_function(
            functions, source=source,
            comment="Batch insert of {} quadratic representatives from {} (source {})".format(
                len(functions), label, source
            )
        )
        print("  inserted entry ids {}..{}".format(ids[0], ids[-1]))
    success("Batch-inserted {} CCZ classes".format(db.number_of_ccz_classes))

```
We finally print a summary of the database, and the journal to check that every source was recorded correctly

```python
    section("Summary")
    print(db)

    section("Journal")
    for entry in db.get_journal():
        print("  [{}] {} — {} — {}".format(
            entry["id"], entry["timestamp"], entry["operation"], entry["comment"]
        ))
```
