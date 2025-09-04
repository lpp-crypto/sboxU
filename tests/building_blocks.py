#!/usr/bin/env sage

from sage.all import *
from sboxUv2 import *
from sage.crypto.sboxes import sboxes


if __name__ == "__main__":
    print("Constructing iScream's Sbox as a Feistel")
    SW=swap_halves(8)
    F1=[0, 8, 6, 13, 5, 15, 7, 12, 4, 14, 2, 3, 9, 1, 11, 10]
    S=feistel_round(F1)
    if SW*(S**3)==Sb(sboxes["iScream"]):
        print("Succeeded")
    else : 
        print("Failed") 