from sboxUv2.core import get_sbox
from sage.all import GF, PolynomialRing

field = GF(2**6)
X = PolynomialRing(field, "X").gen()
g = field.gen()

ccz_class_representatives = [
    # Banff list for quadratic functions
    get_sbox(X**3),
    get_sbox(X**3 + g**11*X**6 + g*X**9),
    get_sbox(g*X**5 + X**9 + g**4*X**17 + g*X**18 + g**4*X**20 + g*X**24 + g**4*X**34 + g*X**40),
    get_sbox(g**7*X**3 + X**5 + g**3*X**9 + g**4*X**10 + X**17 + g**6*X**18),
    get_sbox(X**3 + g*X**24 + X**10), # 4 <-- KIM
    get_sbox(X**3 + g**17*(X**17 + X**18 + X**20 + X**24)),
    get_sbox(X**3 + g**11*X**5 + g**13*X**9 + X**17 + g**11*X**33 + X**48),
    get_sbox(g**25*X**5 + X**9 + g**38*X**12 + g**25*X**18 + g**25*X**36),
    get_sbox(g**40*X**5 + g**10*X**6 + g**62*X**20 + g**35*X**33 + g**15*X**34 + g**29*X**48),
    get_sbox(g**34*X**6 + g**52*X**9 + g**48*X**12 + g**6*X**20 + g**9*X**33 + g**23*X**34 + g**25*X**40),
    get_sbox(X**9 + g**4*(X**10 + X**18 ) + g**9*(X**12 + X**20 + X**40 )),
    get_sbox(g**52*X**3 + g**47*X**5 + g*X**6 + g**9*X**9 + g**44*X**12 + g**47*X**33 + g**10*X**34 + g**33*X**40),
    get_sbox(g*(X**6 + X**10 + X**24 + X**33) + X**9 + g**4*X**17),
    # cubic function
    get_sbox([0, 0, 0, 8, 0, 26, 40, 58, 0, 33, 10, 35, 12, 55, 46, 29, 0, 11, 12, 15, 4, 21, 32, 57, 20, 62, 18, 48, 28, 44, 50, 10, 0, 6, 18, 28, 10, 22, 48, 36, 8, 47, 16, 63, 14, 51, 62, 11, 5, 24, 27, 14, 11, 12, 61, 50, 25, 37, 13, 57, 27, 61, 39, 9])
]


            
