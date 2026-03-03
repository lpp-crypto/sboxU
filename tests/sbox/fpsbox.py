# -*- python -*-


from sboxUv2 import *
from sage.crypto.sboxes import sboxes
from sage.all import *

if __name__ == "__main__":

    with Experiment("Testing S-box module in F_p"):
        for n in range(3, 11):
            print("\n\n----\n", n)
            for t in range(0, 3):
                s = random_function_S_box(randint(2, 7), n)
                b = s.to_bytes()
                print("\n{}\n{}\n{}".format(
                    b,
                    s,
                    Sb(b)))
                
        ### Fp testing session
        p = 3
        Fp = GF(p)
        # Build an SBox from F_3^2 to itself, invertible
        lut = [[0,1],[1,0],[0,2],[0,0],[2,0],[2,2],[1,2],[2,1],[1,1]]
        lut = [[Fp(x),Fp(y)] for x,y in lut]
        u = Sb(lut)
        print(u[(Fp(0),Fp(1))]) 
      
        # R = PolynomialRing(Fp,2,"x")
        # x1, x2 = R.gens()
        # P1 = x1**3 + x2**5 + x1*x2
        # P2 = x2
        # v = Sb([P1,P2])
        # print(u)
        # print(v)
        # a = v.derivative([int(0),int(0)])
        # print(a)
    
        # print(v+u)
        # print(v*u)
        # print(u.inverse())
        # try :
        #     print((v*u).inverse())
        # except Exception as e :
        #     pass
    
        # w = pow(u,0)
        # print(w[(0,1)])
    
