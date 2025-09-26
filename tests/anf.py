from sage.all import *
from sboxUv2 import *
from sage.crypto.sboxes import sboxes
from random import randint


Sb0=Sb(sboxes["Midori_Sb0"])
anf_Sb0=algebraic_normal_form(Sb0)

# the coordinates of the ANF (such as the first one) are SAGE polynomials, meaning their `parent()` is the multivariate polynomial over which they are defined.
R = anf_Sb0[0].parent()
(x0,x1,x2,x3,) = R.gens()

for coord in anf_Sb0:
    print(coord)
pprint(degree_spectrum(Sb0))

if anf_Sb0==[x0*x1*x2 + x0*x1*x3 + x0*x2 + x0*x3 + x1*x2*x3 + x1,
 x0*x2 + x0*x3 + x0 + x2*x3 + x2,
 x0*x1*x2 + x0*x1*x3 + x0*x3 + x0 + x1*x2*x3 + x3 + 1,
 x0*x1*x3 + x0*x1 + x1*x2*x3 + x1*x3 + x2*x3 + 1]:
    print("The ANF of Midori_Sb0 is computed correctly")


test=Sb(anf_Sb0)

if Sb0.lut()==test:
    print("The polynomials can be evaluated to retrieve the LUT")


print("Testing eval_anf on a 8-bit random Boolean Function")
S=Sb([randint(0,1) for _ in range(256)])
anf=algebraic_normal_form(S)[0]
if [eval_anf(anf,x) for x in range(256)]==S:
    print("Success")
else :
    print("Failure")

print("Testing eval_vect_anf on a random 8-bit Sbox")
S=Sb([randint(0,255) for _ in range(256)])
anfs=algebraic_normal_form(S)
if [eval_vect_anf(anfs,x) for x in range(256)]==S:
    print("Success")
else :
    print("Failure")

