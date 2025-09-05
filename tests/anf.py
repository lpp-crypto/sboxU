from sage.all import *
from sboxUv2 import *
from sage.crypto.sboxes import sboxes


Sb0=Sb(sboxes["Midori_Sb0"])
anf_Sb0=algebraic_normal_form(Sb0)
R = PolynomialRing(GF(2), names=('x0', 'x1', 'x2','x3',))
(x0,x1,x2,x3,) = R._first_ngens(4)
for coord in anf_Sb0:
    print(coord)
pprint(degree_spectrum(Sb0))
if anf_Sb0==[x0*x1*x2 + x0*x1*x3 + x0*x2 + x0*x3 + x1*x2*x3 + x1,
 x0*x2 + x0*x3 + x0 + x2*x3 + x2,
 x0*x1*x2 + x0*x1*x3 + x0*x3 + x0 + x1*x2*x3 + x3 + 1,
 x0*x1*x3 + x0*x1 + x1*x2*x3 + x1*x3 + x2*x3 + 1]:
    print("The ANF of Midori_Sb0 is computed correctly")


test=[]
for i in range(16):
    t0,t1,t2,t3=to_bin(i,4)
    test.append(from_bin([anf_Sb0[j](t0,t1,t2,t3) for j in range(4)]))

if Sb0.lut()==test:
    print("The polynomials can be evaluated to retrieve the LUT")
