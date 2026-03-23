from sage.all import *
from sboxUv2 import *


if __name__ == "__main__":
    with Experiment("Testing PRNG"):

        section("Basic queries")
        
        p = InsecurePRNG(bytearray([0, 1]))
        x = 0
        for t in range(0, 2**20):
            x += p(0, 30)
            # x += randint(0, 30)
        print("tot: ", x)

        section("Random S-boxes")

        for i in range(0, 10):
            subsection("new random S-boxes {}".format(i))

            print("bijection")
            
            s = rand_invertible_S_box(p, 5)
            pprint(s)
            pprint(differential_spectrum(s))
            pprint(absolute_walsh_spectrum(s))
            print("")

            print("non-bijection")
            
            s = rand_S_box(p, 5, p(4,7))
            pprint(s)
            pprint(differential_spectrum(s))
            pprint(absolute_walsh_spectrum(s))
            print("")
        
