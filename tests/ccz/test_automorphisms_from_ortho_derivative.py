#!/usr/bin/env python
import sys
from sage.all import *
from sboxU import *
from time import time

DB_PATH = sixBitAPNs()



def main_test():
    with Experiment('automorphisms_from_ortho_derivative — mode benchmark'):
        section('Count and automorphism verification')
        # --- { 
        with APNFunctions(DB_PATH) as db:
            entries = db.query_functions({"degree": 2})
        pprint("Loaded {} quadratic APN functions".format(len(entries)))
        t_std = t_prod = 0.0
        for entry in entries:
            f     = entry["sbox"]
            label = "id={}".format(entry["id"])
            t = time(); std  = automorphisms_from_ortho_derivative(f, mode="standard"); t_std  += time() - t
            t = time(); prod = automorphisms_from_ortho_derivative(f, mode="product");  t_prod += time() - t
            if len(std) != len(prod):
                fail("{}: standard={} automorphisms, product={}".format(label, len(std), len(prod)))
                continue
            for L in prod:
                g = ccz_equivalent_function(f, L)
                if len(g) == 0 or g.lut() != f.lut():
                    fail("{}: map is not a graph automorphism".format(label))
                    break
            else:
                success("{}: {} automorphisms, all verified".format(label, len(prod)))
        pprint("Total time — standard: {:.3f}s  product: {:.3f}s".format(t_std, t_prod))
        # --- } 
    return exit_code()


if __name__ == '__main__':    sys.exit(main_test())
