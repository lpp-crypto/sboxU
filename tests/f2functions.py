#!/usr/bin/env sage

from sage.all import *

from sboxUv2 import *


n = 10

def pretty_bin(x):
    result = ""
    while x != 0:
        result = str(x & 1) + result
        x = x >> 1
    while len(result) <= n:
        result = "0" + result
    return result


if __name__ == "__main__":
    canonical = [1 << i for i in range(0, n)]
    entries = []
    for t in range(0, 100):
        mask = randint(0, 2**n)
        if hamming_weight(mask) < 4:
            x = linear_combination(canonical, mask)
            print("{} {:3d} {:3d}".format(
                pretty_bin(x),
                lsb(x),
                msb(x)))
            entries.append(x)
    print(pretty_bin(xor(entries)) + " (tot)")
    print(pretty_bin(xor(entries, 0x3FF)) + " (not tot)")
