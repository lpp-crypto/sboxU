from sboxUv2.core.sbox import Sb
from sage.all import PolynomialRing,inverse_mod

def closed_butterfly(alpha,beta):
    if alpha.parent() != beta.parent() :
        raise Exception("alpha and beta must belong to the same field")
    else :
        R=PolynomialRing(alpha.parent(), names =('x','y'))
        (x,y,) = R._first_ngens(2)
        poly=(x+alpha*y)**3+beta*(y**3)
        return Sb([poly(x,y),poly(y,x)])
    
def open_butterfly(alpha,beta):
    if alpha.parent() != beta.parent() :
        raise Exception("alpha and beta must belong to the same field")
    else :
        R=PolynomialRing(alpha.parent(), names =('x','y'))
        n=R.base_ring().degree()
        (x,y,) = R._first_ngens(2)
        poly=(x+alpha*y)**3+beta*(y**3)
        poly_inv=(beta*y**3+x)**(inverse_mod(3,2**n-1))+alpha*y
        return Sb([poly(y,poly_inv(x,y)),poly_inv(x,y)])