from sage.all import *
from sboxU import *

from collections import defaultdict
import itertools


def super_walsh(s):
    result = defaultdict(int)
    for b in range(0, len(s)):
        F_a = []
        for x in range(0, len(s)):
            F_a.append( oplus(b, s[x]) )
        w = walsh_spectrum(F_a)
        result[str(w)] += 1
    return result


with Experiment("Testing APN functions-related functions"):

    # section("Initialization")
    # n = 7
    # g = GF(2**n)
    # cube = monomial(3, g)
    # print("n=", n)
    # print(cube)

    # section("Testing mugshot for quadratic APN functions")
    # pprint(differential_spectrum(cube))
    # pprint(ortho_derivative(cube))
    # pprint(sigma_multiplicities(cube, 4))
    # pprint(thickness_spectrum(cube))
    
    # print(apn_ea_mugshot(cube))

    # section("Testing CCZ-class reduction")
    
    # reprs = enumerate_ea_classes_apn_quadratic(cube)
    # print("EA classes in CCZ: ", len(reprs))
    # for s in reprs:
    #     pprint(degree_spectrum(s))
    #     print(s)
    #     pprint("checking: ", degree_spectrum(ccz_equivalent_quadratic_function(s)))

    # section("Testing DB")
    
    with APNFunctions("sboxU/scripts/apnDB/apnDB.db") as db:
        print(db[0])
        s = db[0]["sbox"]
        pprint(s)
        pprint(differential_spectrum(s))
        pprint(thickness_spectrum(s))
        pprint(degree_spectrum(s))
        pprint(absolute_walsh_spectrum(s))
