import sys
from sage.all import *
from sboxU import *
# --- { 
import os, random
from time import time
from sage.all import Matrix, GF, vector
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(SCRIPT_DIR, "../sboxU/sboxU/scripts/apnDB/apn6.db")
with APNFunctions(DB_PATH) as db:
    all_apn = db.query_functions({"degree": 2})
pprint("Loaded {} quadratic APN functions".format(len(all_apn)))
def rand_lin_perm(n):
    while True:
        m = Matrix(GF(2), n, n, [[random.randint(0, 1) for _ in range(n)] for _ in range(n)])
        if m.rank() == n:
            return m
def mat_to_sbox(mat):
    n = mat.ncols()
    lut = [from_bin(mat * vector(to_bin(x, n))) for x in range(1 << n)]
    return get_sbox(lut)
def random_ea_equiv(f):
    """Return A*f*B + C with random linear permutations A, B and linear map C."""
    n = f.get_input_length()
    A = mat_to_sbox(rand_lin_perm(n))
    B = mat_to_sbox(rand_lin_perm(n))
    C = mat_to_sbox(Matrix(GF(2), n, n, [[random.randint(0, 1) for _ in range(n)] for _ in range(n)]))
    return A * f * B + C
random.seed(12)
ccz_reps = {}
for e in all_apn:
    if e["ccz_id"] not in ccz_reps:
        ccz_reps[e["ccz_id"]] = e["sbox"]
chosen = random.sample(sorted(ccz_reps.keys()), 6)
src_funcs = [ccz_reps[i] for i in chosen[:3]]
neg_funcs = [ccz_reps[i] for i in chosen[3:]]
ea_funcs  = [random_ea_equiv(f) for f in src_funcs]
pprint("Source CCZ classes:   {}".format(chosen[:3]))
pprint("Negative CCZ classes: {}".format(chosen[3:]))
# --- } 
def main_test():
    with Experiment(' Correctness'):
        section(' Correctness')
        # --- { 
        n_ok = n_fail = 0
        for i in range(3):
            fi, gi = src_funcs[i], ea_funcs[i]
            ref  = are_ea_equivalent(fi, gi)
            std  = len(ea_mapping_from_vq(fi, gi, mode="standard")) > 0
            prod = len(ea_mapping_from_vq(fi, gi, mode="product"))  > 0
            if ref == std == prod == True:
                success("positive pair {}: ref, standard, product all True".format(i))
                n_ok += 1
            else:
                fail("positive pair {}: ref={} std={} prod={}".format(i, ref, std, prod))
                n_fail += 1
        for i in range(3):
            fi, gi = neg_funcs[i], ea_funcs[i]
            ref  = are_ea_equivalent(fi, gi)
            std  = len(ea_mapping_from_vq(fi, gi, mode="standard")) > 0
            prod = len(ea_mapping_from_vq(fi, gi, mode="product"))  > 0
            if ref == std == prod == False:
                success("negative pair {}: ref, standard, product all False".format(i))
                n_ok += 1
            else:
                fail("negative pair {}: ref={} std={} prod={}".format(i, ref, std, prod))
                n_fail += 1
        # --- } 
        section(' Mode comparison')
        # --- { 
        modes = ["standard", "product"]
        timings = {m: 0.0 for m in modes}
        pairs = list(zip(src_funcs + neg_funcs, ea_funcs + ea_funcs))
        for fi, gi in pairs:
            for m in modes:
                t = time(); ea_mapping_from_vq(fi, gi, mode=m); timings[m] += time() - t
        pprint("Mean time per pair ({} pairs):".format(len(pairs)))
        for m in modes:
            pprint("  {}: {:.3f}s".format(m, timings[m] / len(pairs)))
        # --- } 
    return exit_code()


if __name__ == '__main__':    sys.exit(main_test())
