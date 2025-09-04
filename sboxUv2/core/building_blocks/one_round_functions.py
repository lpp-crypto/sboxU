from sboxUv2.core.f2functions import to_bin, from_bin, oplus # type: ignore
from sboxUv2.core.sbox import Sb # type: ignore

def swap_halves(n):
    ## Return an Sbox on n bits which swaps the most significant half with the least significant one.
    if n%2!=0:
        raise ValueError("n has to be even")
    return Sb([from_bin(to_bin(x,n)[n//2:n]+to_bin(x,n)[0:n//2]) for x in range(2**n)])

def feistel_round(F1):
    ## Return an Sbox S which is the result of a classic round feistel with round function F1
    S1=Sb(F1)
    n=S1.get_input_length()
    res=[]
    for xR in range(2**n):
        for xL in range(2**n):
            res.append(from_bin(to_bin(oplus(xL,S1[xR]),n)+to_bin(xR,n))) 
    return Sb(res)


    
    

