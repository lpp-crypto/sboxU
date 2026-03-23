from sage.all import *
from sboxUv2 import *
from random import randint


print("Testing the ANF part of the Sb factory")
test= True
for _ in range(50): ## Generating random 8-bit Sboxes
    S=get_sbox([randint(0,255) for _ in range(256)])
    if (get_sbox(algebraic_normal_form(S))== S)==False:
        test=False

if test :
    print("Succeeded")
else :
    print("Failed")


## !TODO! Test of the multivariate (non-ANF) part


