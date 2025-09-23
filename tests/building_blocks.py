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

    print("Constructing Scream's Sbox as a Feistel")
    SW=swap_halves(8)
    F1=Sb([0,2,0,0xb,3,0,0,0xa,1,0xe,0,6,0xa,4,5,2])
    F2=Sb([0,2,0xc,7,5,0xf,0xd,6,4,0xe,8,9,3,1,0xb,0xa])
    F3=Sb([2,0,0xb,0,0,3,0xa,0,0xe,1,6,0,4,0xa,2,5])
    if feistel_round(F3)*feistel_round(F2)*feistel_round(F1)*SW==Sb(sboxes["Scream"]):
        print("Succeeded")
    else : 
        print("Failed") 