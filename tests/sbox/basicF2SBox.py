import sys
from sage.all import *
from sboxU import *
def main_test():
    with Experiment(' Basic functionalities of the S_box class'):
        section(' Basic functionalities of the S_box class')
        # --- { 
        s = get_sbox([57,29,43,15,16,12,10,59,37,63,32,23,9,24,31,1,33,
            26,39,38,18,7,58,49,27,54,62,52,4,53,5,25,47,
            11,40,0,50,48,8,42,35,13,3,36,2,44,56,6,28,
            34,17,14,45,30,61,21,51,19,46,60,41,20,55,22])
        # --- } 
        # --- { 
        print("basic list representation:")
        print(s)
        print("\npretty printing the same information:")
        pprint(s)
        # --- } 
        subsection(' Differential properties')
        # --- { 
        a = 1
        D_a_s = derivative(s, 1)
        pprint(D_a_s)
        # --- } 
        # --- { 
        derivative_is_translation_invariant = True
        for x in range(0, 2**s.get_output_length()):
            if D_a_s[x] != D_a_s[oplus(x, 1)]:
                fail("derivative should be identical on x and x+a for all x and a, but it isn't the case for x={}, a={}".format(x, a))
                derivative_is_translation_invariant = False
        if derivative_is_translation_invariant:
            success("sanity check passed: the derivative on a is invariant under translation by a")
        # --- } 
        subsection(' Linear properties')
        subsection(' Boomerang properties')
        subsection(' Polynomial representation')
        section(' Building an S-box')
        subsection(' The case of S-boxes from the literature')
        # --- { 
        s = get_sbox("AES")
        # --- } 
        # --- { 
        for test in [("AES", 4),
                     ("Kuznyechik", 8),
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
        section(' Univariate polynomials')
        section(' Operations on S-boxes')
        subsection(' Composition')
    return exit_code()


if __name__ == '__main__':    sys.exit(main_test())
