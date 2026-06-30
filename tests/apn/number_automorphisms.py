from sage.all import *
from sboxU import *


with Experiment("Testing APN functions-related functions"):
    section("number of automorphisms")

    with APNFunctions("./apnDB.db") as db:
        entries = db.query_functions({"degree" : 2})
        for entry in entries:
            aut = automorphisms_from_ortho_derivative(entry["sbox"])
            l = len(aut)
            if l % 5 == 0:
                subsection("Interesting function")
                s = entry["sbox"]
                print(l, factor(l))
                pprint(absolute_walsh_spectrum(s))
                pprint(thickness_spectrum(s))
                mask = 2**6-1
                for L in aut:
                    interesting = False
                    # for x in range(0, 2**6):
                    #     if (L(x) & mask) != L(x):
                    #         interesting = False
                    #         break
                    # for x in range(0, 2**6):
                    #     if (L(x<<6) & (mask << 6)) != L(x << 6):
                    #         interesting = False
                    #         break
                    if L*L*L*L*L*L == L:
                        interesting = True
                    if interesting:
                        print("\n")
                        print(str(L))
