from sboxUv2.core.f2functions import to_bin, from_bin, oplus
from sboxUv2.core.sbox import get_sbox

def swap_halves(n):
    """
    Args:
       n: The bit length of the result. 

    Returns:
        A n-bit S_box which swaps the most significant half of the input and the least significant one. 
    """
    if n%2!=0:
        raise ValueError("n has to be even")
    return get_sbox([from_bin(to_bin(x,n)[n//2:n]+to_bin(x,n)[0:n//2]) for x in range(2**n)])

def feistel_round(F1):
    """
    Args:
       F1: an S_box-able object. 

    Returns:
        An S_box which is the result of a classic feistel round with round function F1
    """
    S1=get_sbox(F1)
    n=S1.get_input_length()
    res=[]
    for xR in range(2**n):
        for xL in range(2**n):
            res.append(from_bin(to_bin(xR,n)+to_bin(oplus(xL,S1[xR]),n)))
    return get_sbox(res)

    
    

