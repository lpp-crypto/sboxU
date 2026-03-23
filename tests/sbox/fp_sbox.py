# -*- python -*-


from sboxUv2 import *
from sage.crypto.sboxes import sboxes
from sage.all import *

if __name__ == "__main__":
    ### Fp testing session
    p = 3
    Fp = GF(p)
    # Build an SBox from F_3^2 to itself, invertible
    lut = [[0,1],[1,0],[0,2],[0,0],[2,0],[2,2],[1,2],[2,1],[1,1]]
    lut = [[Fp(x),Fp(y)] for x,y in lut]
    u = get_sbox(lut)
    print(u[(Fp(0),Fp(1))]) 
    print(u.coordinate(0))
  
    R = PolynomialRing(Fp,2,"x")
    x1, x2 = R.gens()
    P1 = x1**3 + x2**5 + x1*x2
    P2 = x2
    v = get_sbox([P1,P2])
    print(v)
    print(v.coordinate(0))
    print(u)

    # print(v+u)
    # print(v*u)
    # print(u.inverse())
    # try :
    #     print((v*u).inverse())
    # except Exception as e :
    #     pass

    # w = pow(u,0)
    # print(w[(0,1)])
