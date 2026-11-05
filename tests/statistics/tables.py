#!/usr/bin/env python
import sys
from sage.all import *
from sboxU import *
from collections import defaultdict


def main_test():
    with Experiment("Basic Statistical Properties in the Binary Case"):
        section("Differential properties")
        # --- {
        s = random_permutation_S_box(6)
        # --- }
        subsection("Derivatives")
        # --- {
        D_1_s = derivative(s, 1)
        pprint(D_1_s)
        # --- }
        # --- {
        derivative_is_translation_invariant = True
        for x in range(0, 2 ** s.get_output_length()):
            if D_1_s[x] != D_1_s[oplus(x, 1)]:
                fail(
                    "derivative should be identical on x and x+a for all x and a, but it isn't the case for x={}, a={}".format(
                        x, a
                    )
                )
                derivative_is_translation_invariant = False
        if derivative_is_translation_invariant:
            success(
                "sanity check passed: the derivative on 1 is invariant under translation by 1"
            )
        # --- }
        subsection("DDT")
        # --- {
        d = ddt(s)
        # --- }
        # --- {
        ddt_row = [0 for x in s.input_space()]
        for x in s.input_space():
            ddt_row[D_1_s[x]] += 1
        if ddt_row == d[1]:
            success("The DDT row corresponding to input difference 1 is correct")
        else:
            fail("Problem with the DDT")
        # --- }
        section("Linear properties")
        subsection("Boomerang properties")
        section("Collision properties")
        sprime_lat = [0, 0, 0, 1, 1, 2, 2, 2]
        sprime = get_sbox(sprime_lat)
        pprint(sprime)
        print(f"Collision table of SBox {sprime}:")
        print(collision_table(sprime))
        print(f"Collision spectrum :")
        print(collision_spectrum(sprime))
        section("Testing spectrum equality")
        ssecond_lat = [3, 3, 3, 2, 2, 2, 0, 0]
        ssecond = get_sbox(ssecond_lat)
        # sprime and ssecond should have the same collision spectrum
        coll_spec_sprime = collision_spectrum(sprime)
        coll_spec_ssecond = collision_spectrum(ssecond)
        print(f"First collision spectrum : {coll_spec_sprime}, from S_box : {sprime}")
        print(
            f"Second collision spectrum : {coll_spec_ssecond}, from S_box : {ssecond}"
        )
        print("Does the equality test work on the two spectrums ?")
        if coll_spec_sprime == coll_spec_ssecond:
            success("The equality test works")
        else:
            fail("Problem with the equality test of spectrums")
        section("Testing hashing of spectrums")
        if hash(coll_spec_sprime) == hash(coll_spec_ssecond):
            success("Hashing equal spectrums returns equal hashes")
        else:
            fail("Hashing equal spectrum does not return equal hashes")
        section("References")
    return exit_code()


if __name__ == "__main__":
    sys.exit(main_test())
