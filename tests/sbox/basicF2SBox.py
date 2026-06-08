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
        print("basic list representation")
        print(s)
        print("pretty printing the same information")
        pprint(s)
        # --- } 
        subsection(' Linear properties')
        subsection(' Differential properties')
        subsection(' Boomerang properties')
        subsection(' Polynomial representation')
        section(' Building an S-box')
        subsection(' The case of S-boxes from the literature')
        # --- { 
        s = get_sbox("AES")
        # --- } 
        # --- { 
        for test in [("AES", 4),
                     ("Kuznyechik", 8)]:
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
    return exit_code()


if __name__ == '__main__':    sys.exit(main_test())
