# -*- python -*-


from sboxU import *
from sage.crypto.sboxes import sboxes
from sage.all import *

if __name__ == "__main__":

    with Experiment("Testing Linear Casts in F_2"):

        print("In SboxU, an S-box is seen as a lookup table, i.e., as a function taking integers in \\{0,...,l-1\\} as its input. For functions over F_2, the usual mapping from F_2^n to \\{0,...,2^{n-1}\\} is then implicitely used. However, SboxU allows you to specify other functions to cast elements from a set of your choice to the integers (and back). Then, when using the __call__ method, these casts are used to transform the input and the output of the S-box.")
        
        section("Finite field")

        print("For S-boxes defined from polynomials over finite fields, a cast sending field elements to integer and its inverse are automatically added during the S-box construction. Thus, evaluating the S-box on a field element returns a field element.")
        
        g = GF(2**8)
        X = PolynomialRing(g, "X").gen()
        s = get_sbox((X**16 + X) * X**3)
        for t in range(0, 10):
            x     = g.random_element()
            x_int = ffe_to_int(x)
            print("{:3d} = {:40}  ->  {:40} ({})".format(
                x_int,
                str(x),
                str(s(x)),
                ffe_from_int(s[x_int], g) == s(x)
            ))


        section("Cartesian Products")

        subsection("Principle")
        
        print("The input of an S-box is sometimes better described as an element of F_2^t×F_2^{n-t}, and similarily for its output. To better capture this, it is possible to define casts sending cartesian products to integers (and vice versa).")

        n = 5
        t = 2
        s = random_permutation_S_box(5)
        pprint(s)


        c     = CastFromF2Product([t, n-t])
        c_inv = CastToF2Product([t, n-t])
        s.attach_casts_pair(c, c_inv)
        line = "  "
        for x1 in range (0, 2**t):
            line += "    {}   ".format(x1)
        print(line)
        for x2 in range (0, 2**(n-t)):
            line = "{:x} | ".format(x2)
            for x1 in range (0, 2**t):
                line += "{}  ".format(s((x1, x2)))
            print(line)
            
        subsection("Querying multiple values")

        print("It is possible to \"saturate\" an input and to get the correponding list of inputs using the `structure` method of CastFromF2Product")

        for x in c.structure([2, "*"]):
            print(x, s(x))
