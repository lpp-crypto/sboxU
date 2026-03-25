from sboxU.core.sbox import get_sbox
from sage.all import PolynomialRing,inverse_mod

def closed_butterfly(alpha, beta):
    """Returns the closed butterfly S-box defined by the parameters alpha and beta.

    The closed butterfly is the vectorial Boolean function (x, y) -> (f(x,y), f(y,x))
    where f(x, y) = (x + alpha*y)^3 + beta*y^3, operating on pairs of elements
    from the same finite field of characteristic 2.

    It was introduced in [C:PerUdoBir16], and then refined in [IEEE:CanDuvPer17] (where
    conditions on alpha and beta for APN-ness were identified).

    Args:
        alpha: A finite field element.
        beta: A finite field element belonging to the same field as alpha.

    Returns:
        An S_box instance whose LUT encodes the closed butterfly.

    Raises:
        Exception: If alpha and beta do not belong to the same field.
    """
    if alpha.parent() != beta.parent() :
        raise Exception("alpha and beta must belong to the same field")
    else :
        R=PolynomialRing(alpha.parent(), names =('x','y'))
        (x,y,) = R._first_ngens(2)
        poly=(x+alpha*y)**3+beta*(y**3)
        return get_sbox([poly(x,y),poly(y,x)])
    
def open_butterfly(alpha, beta):
    """Returns the open butterfly S-box defined by the parameters alpha and beta.

    The open butterfly is the vectorial Boolean function (x, y) -> (f(y, g(x,y)), g(x,y))
    where f(x, y) = (x + alpha*y)^3 + beta*y^3 and g(x, y) = (beta*y^3 + x)^{1/3} + alpha*y,
    operating on pairs of elements from the same finite field of characteristic 2 with odd degree.

    It was introduced in [C:PerUdoBir16], and then refined in [IEEE:CanDuvPer17] (where
    conditions on alpha and beta for APN-ness were identified).

    Args:
        alpha: A finite field element.
        beta: A finite field element belonging to the same field as alpha.

    Returns:
        An S_box instance whose LUT encodes the open butterfly.

    Raises:
        Exception: If alpha and beta do not belong to the same field.
    """
    if alpha.parent() != beta.parent() :
        raise Exception("alpha and beta must belong to the same field")
    else :
        R=PolynomialRing(alpha.parent(), names =('x','y'))
        n=R.base_ring().degree()
        (x,y,) = R._first_ngens(2)
        poly=(x+alpha*y)**3+beta*(y**3)
        poly_inv=(beta*y**3+x)**(inverse_mod(3,2**n-1))+alpha*y
        return get_sbox([poly(y,poly_inv(x,y)),poly_inv(x,y)])
