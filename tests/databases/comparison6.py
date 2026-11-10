"""Compare the compact 6-bit APN database against the reference apn6.db.

Checks:
  1. Every quadratic CCZ class (ccz_id 0–12) has the same number of EA
     representatives in both databases.
  2. For CCZ class 6, the degree distribution is identical across both databases.
  3. For every function in CCZ class 6 of the compact database, there is exactly
     one EA-equivalent function in the same class of the reference database.
     Mugshot pre-filtering is used to limit the number of EA-equivalence calls.

Expected runtime: a few seconds (for the per-function equivalence checks in step 3).

Prerequisites: Having compact6.db from database_tutorial_6_bit_compact.py

"""

import os
import sys

from sboxU import *
from sboxU.display import Experiment, section, subsection, success, fail, exit_code
from sboxU.apn import APNFunctions, APNFunctions_compact, sixBitAPNs
from sboxU.apn import apn_ea_mugshot, apn_ea_mugshot_from_spectra, sigma_multiplicities
from sboxU.core import algebraic_degree, degree_spectrum
from sboxU.statistics import absolute_walsh_spectrum
from sboxU.ccz import are_ea_equivalent, thickness_spectrum, are_ea_equivalent_from_quadratic_representative

DB_COMPACT  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "compact6.db")
DB_CLASSIC  = sixBitAPNs()
CCZ_CLASS_FOR_DETAILED_CHECK = 6


def raw_mugshot(sb):
    """Compute the raw (non-truncated) EA mugshot for a given S-box."""
    if algebraic_degree(sb) == 2:
        return bytes(apn_ea_mugshot(sb))
    abs_walsh = absolute_walsh_spectrum(sb)
    deg_spec  = degree_spectrum(sb)
    thk_spec  = thickness_spectrum(sb)
    return bytes(apn_ea_mugshot_from_spectra(
        abs_walsh, deg_spec, sigma_multiplicities(sb, k=4), thk_spec
    ))


def compare():
    with Experiment("Comparing compact DB against reference apn6.db"):

        section("Loading databases")
        print("Compact DB : {}".format(DB_COMPACT))
        print("Classic DB : {}".format(DB_CLASSIC))

        with APNFunctions_compact(DB_COMPACT, n=6) as compact_db, \
             APNFunctions(DB_CLASSIC) as classic_db:

            n_ccz = compact_db.number_of_ccz_classes
            print("Compact — {} entries, {} CCZ classes".format(len(compact_db), n_ccz))
            print("Classic — {} entries, {} CCZ classes".format(
                len(classic_db), classic_db.number_of_ccz_classes
            ))

            # ---------------------------------------------------------------- #
            section("Step 1 — CCZ-class sizes")
            # ---------------------------------------------------------------- #

            size_errors = []
            for ccz_id in range(n_ccz):
                n_compact = len(compact_db.query_functions({"ccz_id": ccz_id}))
                n_classic  = len(classic_db.query_functions({"ccz_id": ccz_id}))
                status = "OK" if n_compact == n_classic else "MISMATCH"
                print("  ccz_id {:2d}: compact={:3d}  classic={:3d}  {}".format(
                    ccz_id, n_compact, n_classic, status
                ))
                if n_compact != n_classic:
                    size_errors.append(ccz_id)

            if size_errors:
                fail("Size mismatch in CCZ classes: {}".format(size_errors))
            else:
                success("All {} CCZ classes have matching sizes".format(n_ccz))

            # ---------------------------------------------------------------- #
            section("Step 2 — Degree distribution in CCZ class {}".format(
                CCZ_CLASS_FOR_DETAILED_CHECK
            ))
            # ---------------------------------------------------------------- #

            compact_class = compact_db.query_functions({"ccz_id": CCZ_CLASS_FOR_DETAILED_CHECK})
            classic_class  = classic_db.query_functions({"ccz_id": CCZ_CLASS_FOR_DETAILED_CHECK})

            compact_degrees = sorted(e["degree"] for e in compact_class)
            classic_degrees  = sorted(e["degree"] for e in classic_class)

            print("  Compact degrees : {}".format(compact_degrees))
            print("  Classic degrees : {}".format(classic_degrees))

            if compact_degrees == classic_degrees:
                success("Degree distributions match")
            else:
                fail("Degree distributions differ")

            # ---------------------------------------------------------------- #
            section("Step 3 — Per-function EA equivalence in CCZ class {}".format(
                CCZ_CLASS_FOR_DETAILED_CHECK
            ))
            print("({} functions in compact class, using mugshot pre-filter)".format(
                len(compact_class)
            ))
            # ---------------------------------------------------------------- #

            subsection("Building mugshot index for classic class")
            classic_by_mugshot = {}
            for entry in classic_class:
                mug = raw_mugshot(entry["sbox"])
                classic_by_mugshot.setdefault(mug, []).append(entry)
            print("  {} distinct mugshots indexed".format(len(classic_by_mugshot)))

            subsection("Checking each compact-DB function against classic DB")
            equiv_errors = []
            for i, c_entry in enumerate(compact_class):
                c_sb       = c_entry["sbox"]
                mug        = raw_mugshot(c_sb)
                candidates = classic_by_mugshot.get(mug, [])
                ea_check = lambda f, g: are_ea_equivalent_from_quadratic_representative(f.lut(), g.lut())
    
                matches    = sum(
                    1 for cl_entry in candidates
                    if ea_check(c_sb, cl_entry["sbox"])
                )
                status = "OK" if matches == 1 else "FAIL (matches={})".format(matches)
                print("  [{:3d}/{}] compact id={:4d}  degree={}  candidates={}  {}".format(
                    i + 1, len(compact_class),
                    c_entry["id"], c_entry["degree"], len(candidates), status
                ))
                if matches != 1:
                    equiv_errors.append(c_entry["id"])

            if equiv_errors:
                fail("EA-equivalence check failed for compact entry ids: {}".format(
                    equiv_errors
                ))
            else:
                success("Every function in CCZ class {} has exactly one "
                        "EA-equivalent in the classic DB".format(
                    CCZ_CLASS_FOR_DETAILED_CHECK
                ))

    sys.exit(exit_code())


if __name__ == "__main__":
    compare()
