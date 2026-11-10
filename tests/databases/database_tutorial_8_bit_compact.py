#!/usr/bin/env python
import sys
from sage.all import *
from sboxU import *
import os

import sys

from sboxU.apn import APNFunctions_compact

from sboxU.scripts.apnDB.reprs8 import (

    all_WenTanGon, first_QAMs, second_QAMs, all_BeiLea, all_BLP22,

)

# (function group getter, human label, biblio_8 source index)

GROUPS = [

    (all_WenTanGon, "Weng, Tan & Gong (2013)",           5),

    (first_QAMs,    "Yu, Wang & Li (2014) - first QAM",  1),

    (second_QAMs,   "Yu & Perrin (2022) - second QAM",   2),

    (all_BeiLea,    "Beierle & Leander (2022)",           3),

    (all_BLP22,     "Beierle, Leander & Perrin (2022)",   4),

]



def main_test():
    with Experiment('Tutorial - Generating the first compact 8-bit APN database'):
        section('Setup')
        # --- { 
        DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "compact8.db")
        if os.path.exists(DB_PATH):
            print("Removing existing database: {}".format(DB_PATH))
            os.remove(DB_PATH)
        print("Database will be written to: {}".format(DB_PATH))
        # --- } 
        section('Insert database')
        # --- { 
        with APNFunctions_compact(DB_PATH, n=8) as db:
        # --- } 
        # --- { 
            section("Batch-inserting quadratic CCZ-class representatives")
            for group_fn, label, source in GROUPS:
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
        # --- } 
        # --- { 
            section("Summary")
            print(db)
            section("Journal")
            for entry in db.get_journal():
                print("  [{}] {} — {} — {}".format(
                    entry["id"], entry["timestamp"], entry["operation"], entry["comment"]
                ))
        # --- } 
    return exit_code()


if __name__ == '__main__':    sys.exit(main_test())
