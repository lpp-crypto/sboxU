from sage.all import *
from sboxUv2 import *


if __name__ == "__main__":
    with Experiment("Testing PRNG"):

        section("Basic queries")
        
        p = UnsafePRNG(bytearray([0, 1]))
        x = 0
        for t in range(0, 2**20):
            x += p(0, 30)
            # x += randint(0, 30)
        print("tot: ", x)

        section("Random S-boxes")

        for i in range(0, 10):
            subsection("new random S-box {}".format(i))
            s = rand_invertible_S_box(p, 5)
            pprint(s)
            pprint(differential_spectrum(s))
            pprint(absolute_walsh_spectrum(s))
            print("")
        
