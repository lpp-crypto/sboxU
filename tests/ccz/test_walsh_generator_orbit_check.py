import sys
from sage.all import *
from sboxU import *
# --- { 
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH    = os.path.join(SCRIPT_DIR, "../sboxU/sboxU/scripts/apnDB/apn6.db")
def affine_key(L):
    return tuple(L.get_S_box().lut())
def group_closure(n, gens):
    """BFS closure of <gens> inside the group of (2n)x(2n) F2AffineMap, via right-multiplication."""
    e = identity_F2AffineMap(2*n)
    elements = {affine_key(e): e}
    frontier = [e]
    while frontier:
        new_frontier = []
        for g in frontier:
            for h in gens:
                gh = g * h
                k = affine_key(gh)
                if k not in elements:
                    elements[k] = gh
                    new_frontier.append(gh)
        frontier = new_frontier
    return elements
# --- } 
def main_test():
    with Experiment(' Derivative automorphisms: canonical-basis generators vs full group'):
        section(' Derivative automorphisms: canonical-basis generators vs full group')
        # --- { 
        with APNFunctions(DB_PATH) as db:
            entries = db.query_functions({"degree": 2})
        pprint("Loaded {} quadratic APN functions".format(len(entries)))
        n_ok = 0
        for entry in entries:
            f = entry["sbox"]
            full_keys = set(affine_key(g) for g in graph_automorphisms_from_derivatives(f))
            gens = gen_set_graph_automorphisms_from_derivatives(f)
            ok = set(group_closure(f.get_input_length(), gens).keys()) == full_keys
            (success if ok else fail)("id={}: {} generators, group size {}, ok={}".format(entry["id"], len(gens), len(full_keys), ok))
            n_ok += ok
        if n_ok == len(entries):
            success("All {} functions: canonical-basis generators close into the full derivative automorphism group".format(n_ok))
        # --- } 
        section(' gen_set_F2AffineMap_group (deterministic) on EL automorphism groups')
        # --- { 
        n_ok = 0
        for entry in entries:
            f = entry["sbox"]
            full_group = graph_el_automorphisms_from_ortho_derivative(f)
            full_keys = set(affine_key(g) for g in full_group)
            gens = gen_set_F2AffineMap_group(full_group, mode="deterministic")
            ok = set(group_closure(f.get_input_length(), gens).keys()) == full_keys
            (success if ok else fail)("id={}: {} generators (deterministic), EL group size {}, ok={}".format(entry["id"], len(gens), len(full_keys), ok))
            n_ok += ok
        if n_ok == len(entries):
            success("All {} functions: gen_set_F2AffineMap_group(mode='deterministic') generates the EL automorphism group".format(n_ok))
        # --- } 
        section(' gen_set_F2AffineMap_group (probabilistic) success rate, CCZ class 0')
        # --- { 
        class0 = next(e for e in entries if e["ccz_id"] == 0)
        f0 = class0["sbox"]
        full_group0 = graph_el_automorphisms_from_ortho_derivative(f0)
        full_keys0 = set(affine_key(g) for g in full_group0)
        pprint("CCZ class 0: |Aut_EL| = {}".format(len(full_group0)))
        n_trials = 10
        n_success = sum(set(group_closure(f0.get_input_length(), gen_set_F2AffineMap_group(full_group0, mode="probabilistic")).keys()) == full_keys0 for _ in range(n_trials))
        rate = n_success / n_trials
        pprint("Probabilistic mode generated the full group in {}/{} trials ({:.1%})".format(n_success, n_trials, rate))
        (success if rate >= 0.5 else fail)("Probabilistic mode success rate: {:.1%} over {} trials".format(rate, n_trials))
        # --- } 
        section(' Walsh_zero_orbits: orbit count vs EA-class count')
        # --- { 
        n_ok = 0
        for entry in entries:
            f = entry["sbox"]
            ws = get_WalshZeroesSpaces(f)
            gens = [B.transpose() for B in automorphisms_from_ortho_derivative(f)]
            orbits = ws.Walsh_zero_orbits(gens)
            n_classes = len(enumerate_ea_classes_apn_quadratic(f))
            partition_ok = sorted(i for o in orbits for i in o) == list(range(len(ws.get_bases())))
            ok = len(orbits) == n_classes and partition_ok
            (success if ok else fail)("id={}: {} orbits, {} EA classes, partition_ok={}".format(entry["id"], len(orbits), n_classes, partition_ok))
            n_ok += ok
        if n_ok == len(entries):
            success("All {} functions: Walsh_zero_orbits count matches EA-class count, valid partition".format(n_ok))
        # --- } 
        section(' enumerate_ea_classes_apn_quadratic: standard vs product vs generators vs product_generator')
        # --- { 
        from time import time
        
        modes = ["standard", "product", "generators"]
        t = {m: 0.0 for m in modes}
        n_ok = 0
        for entry in entries:
            f = entry["sbox"]
            counts = {}
            for m in modes:
                t0 = time(); counts[m] = len(enumerate_ea_classes_apn_quadratic(f, mode=m)); t[m] += time() - t0
            ok = len(set(counts.values())) == 1
            (success if ok else fail)("id={}: {}".format(entry["id"], counts))
            n_ok += ok
        if n_ok == len(entries):
            success("All {} functions: standard, product, generators agree on EA-class count".format(n_ok))
        for m in modes:
            pprint("  {}: {:.3f}s".format(m, t[m]))
                # --- } 
    return exit_code()


if __name__ == '__main__':    sys.exit(main_test())
