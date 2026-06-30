#!/usr/bin/env python
import sys
from sage.all import *
from sboxU import *


def main_test():
    with Experiment('Dealing with S-boxes in F_2'):
        section('Basic functionalities of the S_box class')
        # --- { 
        N = 6
        s = get_sbox([5*x % 2**N for x in range(0, 2**N)])
        # --- } 
        # --- { 
        print("basic list representation:")
        print(s)
        print("\npretty printing the same information:")
        pprint(s)
        # --- } 
        subsection('Queries')
        # --- { 
        if s[0] == 0:
            success("the output of s on the all zero bit string is indeed {}".format(s[0]))
        else:
            fail("s[0] should be 0, something went wrong")
        # --- } 
        # --- { 
        if s[1] == 5:
            success("the output of s on (1,0,...,0) is indeed {}".format(s[0]))
        else:
            fail("s[1] should be 5, something went wrong")
        # --- } 
        subsection('Polynomial representation')
        section('Building an S-box')
        subsection('The case of S-boxes from the literature')
        # --- { 
        s = get_sbox("AES")
        # --- } 
        # --- { 
        for test in [("AES", 4),
                     ("PRESENT", 4),
                     ("Ascon", 8)]:
            name, expected_uniformity = test
            u = differential_uniformity(get_sbox(name))
            if u == expected_uniformity:
                success("{} S-box differential uniformity is indeed {}".format(
                    name,
                    expected_uniformity
                ))
            else:
                fail("{} S-box differential uniformity should be {}, not {}".format(
                    name,
                    expected_uniformity,
                    u
                ))
        # --- } 
        section('Univariate polynomials')
        # --- { 
        field = GF(2**5)
        g = field.gen()
        X = PolynomialRing(field, "x").gen()
        cube = get_sbox(X**3)
        cube_plus = get_sbox(X**3 + g*X)
        # --- } 
        section('Operations on S-boxes')
        subsection('Addition')
        # --- { 
        diff = cube + cube_plus 
        is_linear = (differential_uniformity(diff) == 2**5)
        if is_linear:
            success("the sum of X^3 and X^3+gX is indeed an affine function")
        else:
            fail("something went wrong, gX should be affine")
        # --- } 
        subsection('Composition')
        section('References')
    return exit_code()


if __name__ == '__main__':    sys.exit(main_test())
