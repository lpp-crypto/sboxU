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

    print("Constructing ZUC_S0 Sbox as a Feistel")
    F1=Sb([9,15,0,14,15,15,2,10,0,4,0,12,7,5,3,9])
    F2=Sb([8,13,6,5,7,0,12,4,11,1,14,10,15,3,9,2])
    F3=Sb([2,6,10,6,0,13,10,15,3,3,13,5,0,9,12,13])
    if circ_shift_F2AffineMap(8,-5).get_S_box()*feistel_round(F3)*feistel_round(F2)*feistel_round(F1)*SW==Sb(sboxes["ZUC_S0"]):
        print("Succeeded")
    else :
        print("Failed")


    # print("Testing the construction of closed generalised butterflies")

    # ## Testing Main Theorem of [CDP17]
    # test= True

    # for alpha in GF(8).list()[1::]:
    #     for beta in GF(8).list()[1::]:
    #         if algebraic_degree(closed_butterfly(alpha,beta))!=2:
    #             test=False 

    # for alpha in GF(8).list()[1::]:
    #     if alpha.trace()==0:
    #         for beta in [alpha**3+alpha,alpha**3+alpha.inverse()]:
    #             butterfly=closed_butterfly(alpha,beta)
    #             if differential_uniformity(butterfly) > 2 or linearity(butterfly) != 16:
    #                 test=False

    # for n in [3,5]:     
    #     for alpha in GF(2**n).list()[1::]:
    #         beta=(1+alpha)**3
    #         if beta !=0:
    #             butterfly=closed_butterfly(alpha,beta)
    #             if differential_uniformity(butterfly)!= 2**(n+1) or linearity(butterfly) != 2**((3*n+1)//2):
    #                 test=False 
    # if test :
    #     print("Succeeded")
    # else :
    #     print("Failed")

    # print("Checking if open generalised butterflies are equal to Feistel's")
   
    # for n in [3,5]:
    #     test=True
    #     print("Testing n=",n)
    #     SW=swap_halves(2*n)
    #     R = PolynomialRing(GF(2**n), names=('X')); 
    #     (X,) = R._first_ngens(1)
    #     elts=GF(2**n).list()[1::]

    #     for beta in elts:
    #         if (open_butterfly(GF(2**n)(1),beta) == SW *feistel_round(Sb(beta*X**3))*feistel_round(Sb(X**inverse_mod(3,2**n-1)))*feistel_round(Sb(beta*X**3)))==False:
    #             test=False
    #     if test :
    #         print("Succeeded")
    #     else :
    #         print("Failed")

    print("Constructing Midori128 Sboxes")
    Sb1=Sb(sboxes["Midori_Sb1"])
    permutations=[bit_permutation_F2AffineMap(p).get_S_box() for p in [[4,1,6,3,0,5,2,7],[1,6,7,0,5,2,3,4],[2,3,4,1,6,7,0,5],[7,4,1,2,3,0,5,6]]]
    SSb=[p.inverse()*(Sb1 | Sb1)*p for p in permutations]
    if ([(S.inverse(), linearity(S),differential_uniformity(S)) == (S, 128,64) for S in SSb] ==[True for _ in range(4)]):
        print("Succeeded")
    else: 
        print("Failed")