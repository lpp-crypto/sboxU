from sboxUv2.core.f2functions import get_F2AffineMap,rank_of_vector_set
from sboxUv2.core.sbox import  get_sbox
from random import randint,shuffle



def rand_linear_permutation(n):
    """Returns a random F2AffineMap which is a permutation of {0,1}^N. This is done by drawing vectors at random and adding them to the image if they are not in the span of the current image. 
    """
    borne_max= 2**n -1
    basis=[]
    r = 0
    while r< n:
        x=randint(1,borne_max)
        new_basis= basis + [x]
        new_r = rank_of_vector_set(new_basis)
        if new_r > r :
            basis =new_basis[:]
            r =new_r
    return get_F2AffineMap(basis)


def rand_linear_function(n,m):
    """Returns a random F2AffineMap from {0,1}^n to {0,1}^m. This is done by drawing the vectors of the image at random."""
    image=[]
    borne_max= 2**m-1
    for _ in range(n):
        image.append(randint(0,borne_max))
    return get_F2AffineMap(image)

def rand_Sbox(n,m):
    """Returns a random n bits to m bits S_box"""
    borne_max= 2**m-1
    return get_sbox([randint(0,borne_max) for _ in range(2**n)])

def rand_invertible_Sbox(n):
    """Returns a permutation of {0,1}^n"""
    res=[i for i in range(2**n)]
    shuffle(res)
    return get_sbox(res)



