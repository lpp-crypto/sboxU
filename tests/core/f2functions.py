#!/usr/bin/env python
import sys
from sage.all import *
from sboxU import *
def pretty_bin(x, n):

    return "".join(str(b) for b in reversed(to_bin(x, n)))

    

prng = InsecurePRNG(b"seed")



def main_test():
    with Experiment('Basic Operations on F_2 Vectors'):
        section('Converting between integers and bit lists: `to_bin` and `from_bin`')
        # --- { 
        n = 8
        x = 0b01101001
        print("x           =", pretty_bin(x, n))
        print("to_bin(x,n) =", to_bin(x, n))
        # --- } 
        # --- { 
        from sage.all import randint
        ok = True
        for i in range(1000):
            y = prng(0, 2**n)
            if from_bin(to_bin(y, n)) != y:
                ok = False
        if not ok:
            fail("to_bin/from_bin round-trip")
        # --- } 
        section('Extremal bits: `lsb` and `msb`')
        # --- { 
        for x in [0b00000001, 0b00001000, 0b01101001, 0]:
            print("x = {:>8s}   lsb(x) = {}   msb(x) = {}".format(
                pretty_bin(x, n), lsb(x), msb(x)
            ))
        # --- } 
        # --- { 
        for _ in range(10):
            y = randint(1, 2**n - 1)
            assert (y >> msb(y)) == 1
        print("x >> msb(x) == 1 for all tested non-zero x")
        # --- } 
        section('Counting set bits: `hamming_weight`')
        # --- { 
        for x in [0, 1, 0b01101001, 2**n - 1]:
            print("x = {:>8s}   hamming_weight(x) = {}".format(pretty_bin(x, n), hamming_weight(x)))
        # --- } 
        # --- { 
        for _ in range(10):
            y = randint(0, 2**n - 1)
            assert hamming_weight(y) == sum(to_bin(y, n))
        print("hamming_weight(x) == sum(to_bin(x,n)) for all tested x")
        # --- } 
        section('The canonical scalar product: `scal_prod`')
        # --- { 
        x = 0b01101001
        y = 0b01000011
        print("x            =", pretty_bin(x, n))
        print("y            =", pretty_bin(y, n))
        print("x AND y      =", pretty_bin(x & y, n))
        print("scal_prod(x,y) =", scal_prod(x, y))
        # --- } 
        # --- { 
        for i in range(10):
            a = prng(0, 2**n)
            b = prng(0, 2**n)
            assert scal_prod(a, b) == hamming_weight(a & b) % 2
        print("scal_prod(x,y) == hamming_weight(x & y) % 2 for all tested x,y")
        # --- } 
        section('Comments')
    return exit_code()


if __name__ == '__main__':    sys.exit(main_test())
