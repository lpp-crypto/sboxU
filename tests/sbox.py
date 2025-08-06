# -*- python -*-

from sboxUv2 import *
from sage.crypto.sboxes import sboxes
from sage.all import GF, Integer

if __name__ == "__main__":
    u = Sb(list(range(0, 16)))
    s = random_permutation_S_box(4)
    t = random_function_S_box(4, 2, name="t")
    s_prime = Sb(list(s), name="S'")
    print("u      |", u)
    print("t      |", t)
    print("s      |", s)
    for i in range(0, s.get_input_length()):
        print(s.coordinate(i))
        print(s.component(1 << i))

    print("lut(s) |", s.lut())
    print("t == u |", t == u)
    print("t == s |", t == s)
    print("s' == s|", s_prime == s)
    
    print("s[2] = {}, t[1] = {}, u[2] = {}".format(s[2], t[1], u[2]))
    print("input space for u", list(u.input_space()))

    print("hash(s) =       |", hash(s))
    print("t == [0,1,2,3]  |", t == [0,1,2,3])
    print("t == list(t)    |", t == list(t))


    print("s + t           |", s + t)
    print("s is invertible |", s.is_invertible())
    print("inverse of s: ", s.inverse())
    sq = s**2
    sq.rename("S**2")
    print("s**2", sq)
    print("s**-1: ", s**-1)
    print("s**-1 == s.inverse(): ", s.inverse() == s**-1)
    print("s+t is invertible", (s+t).is_invertible())
    print("can we invert s+t?")
    try:
        print((s+t).inverse())
    except:
        print("no, it indeed threw an exception")

    print("t o s", t * s)

    mul = F2_mul(2, GF(16))
    print(mul)
    comp = s * mul
    print("comp = s o mul", comp)
    print("mul**-1 o s**-1", (mul**-1) * (s**-1))
    print("comp**-1", comp.inverse())
    print("(mul**-1 o s**-1)**-1", ((mul**-1) * (s**-1))**-1)

    print("+4", F2_trans(0x3, bit_length=4))

    g = GF(16)
    print("gf_inv", monomial(14, g))

    sb0 = Sb(sboxes["Midori_Sb0"], name="Midori_Sb0")
    print("Midori_Sb0", sb0)
    print("Midori_Sb0 ** 2", sb0 ** 2)
    
    for delta in sb0.input_space():
        print(sb0.derivative(delta))
