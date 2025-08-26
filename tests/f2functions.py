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
    for t in range(0, 10):
        mask = randint(0, 2**n-1)
        if hamming_weight(mask) < 4:
            x = linear_combination(canonical, mask)
            print("{} {:3d} {:3d}".format(
                pretty_bin(x),
                lsb(x),
                msb(x)))
            entries.append(x)
    print(pretty_bin(xor(entries)) + " (tot)")
    print(pretty_bin(xor(entries, 0x3FF)) + " (not tot)")
            
    tot_sum  = zero_BinLinearMap(n)
    tot_prod = identity_BinLinearMap(n)
    for t in range(0, 10):
        masks = [randint(1, 2**n-1) for i in range(0, randint(3, n+1))]
        print(masks)
        L = Blm(masks)
        x = [randint(0, 2**n) for u in range(0, 15)]
        img = [L(x_i) for x_i in x]
        print(img)
        print([linear_combination(masks, x_i) for x_i in x])
        print(rank_of_vector_set(img), L.rank())
        print("")
        if L.get_input_length() == n:
            tot_prod = L * tot_prod
            tot_sum = L + tot_sum
    print("\ntot_prod", tot_prod)
    print("\ntot_sum", tot_sum)

    ## Testing to_bin and from_bin
    is_correct=True
    for _ in range(10):
        x=randint(0,2**64-1)
        if from_bin(to_bin(x,64))!=x:
            is_correct=False
    if is_correct:
        print("to_bin and from_bin seem to work")
    else :
        print("Error in to_bin or in from_bin")