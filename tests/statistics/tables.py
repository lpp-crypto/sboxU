#!/usr/bin/env python
import sys
from sage.all import *
from sboxU import *
from collections import defaultdict



def main_test():
    with Experiment('Basic Statistical Properties in the Binary Case'):
        section('Differential properties')
        # --- { 
        s = random_permutation_S_box(6)
        # --- } 
        subsection('Derivatives')
        # --- { 
        D_1_s = derivative(s, 1)
        pprint(D_1_s)
        # --- } 
        # --- { 
        derivative_is_translation_invariant = True
        for x in range(0, 2**s.get_output_length()):
            if D_1_s[x] != D_1_s[oplus(x, 1)]:
                fail("derivative should be identical on x and x+a for all x and a, but it isn't the case for x={}, a={}".format(x, a))
                derivative_is_translation_invariant = False
        if derivative_is_translation_invariant:
            success("sanity check passed: the derivative on 1 is invariant under translation by 1")
        # --- } 
        subsection('DDT')
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
        section('Linear properties')
        subsection('Boomerang properties')
        section('References')
    return exit_code()


if __name__ == '__main__':    sys.exit(main_test())
